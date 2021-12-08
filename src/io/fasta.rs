// Copyright 2014-2018 Johannes Köster, Christopher Schröder, Henning Timm.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.
//! Structs and trait to read and write files in FASTA format.
//!
//! # Example
//!
//! ## Read
//!
//! In this example, we parse a fasta file from stdin and compute some statistics
//!
//! ```
//! use bio::io::fasta;
//! use std::io;
//!
//! let mut reader = fasta::Reader::new(io::stdin());
//!
//! let mut nb_reads = 0;
//! let mut nb_bases = 0;
//!
//! for result in reader.records() {
//!     let record = result.expect("Error during fasta record parsing");
//!     println!("{}", record.id());
//!
//!     nb_reads += 1;
//!     nb_bases += record.seq().len();
//! }
//!
//! println!("Number of reads: {}", nb_reads);
//! println!("Number of bases: {}", nb_bases);
//! ```
//!
//! We can also use a `while` loop to iterate over records.
//! This is slightly faster than the `for` loop.
//! ```
//! use bio::io::fasta;
//! use std::io;
//! let mut records = fasta::Reader::new(io::stdin()).records();
//!
//! let mut nb_reads = 0;
//! let mut nb_bases = 0;
//!
//! while let Some(Ok(record)) = records.next() {
//!     nb_reads += 1;
//!     nb_bases += record.seq().len();
//! }
//!
//! println!("Number of reads: {}", nb_reads);
//! println!("Number of bases: {}", nb_bases);
//! ```
//!
//! ## Write
//!
//! In this example we generate 10 random sequences with length 100 and write them to stdout.
//!
//! ```
//! use std::io;
//! use bio::io::fasta;
//!
//! let mut seed = 42;
//!
//! let nucleotides = [b'A', b'C', b'G', b'T'];
//!
//! let mut writer = fasta::Writer::new(io::stdout());
//!
//! for _ in 0..10 {
//!     let seq = (0..100).map(|_| {
//!         seed = ((seed ^ seed << 13) ^ seed >> 7) ^ seed << 17; // don't use this random generator
//!         nucleotides[seed % 4]
//!     }).collect::<Vec<u8>>();
//!
//!    writer.write("random", None, seq.as_slice()).expect("Error writing record.");
//! }
//! ```
//!
//! ## Read and Write
//!
//! In this example we filter reads from stdin on sequence length and write them to stdout
//!
//! ```
//! use bio::io::fasta;
//! use bio::io::fasta::FastaRead;
//! use std::io;
//!
//! let mut reader = fasta::Reader::new(io::stdin());
//! let mut writer = fasta::Writer::new(io::stdout());
//! let mut record = fasta::Record::new();
//!
//! while let Ok(()) = reader.read(&mut record) {
//!     if record.is_empty() {
//!         break;
//!     }
//!
//!     if record.seq().len() > 100 {
//!         writer
//!             .write_record(&record)
//!             .ok()
//!             .expect("Error writing record.");
//!     }
//! }
//! ```


use std::cmp::min;
use std::collections;
use std::convert::AsRef;
use std::fs;
use std::io;
use std::io::prelude::*;
use std::path::Path;
use serde::{Serialize, Deserialize};

use anyhow::Context;
use std::fmt;

/// Type alias for an owned text, i.e. ``Vec<u8>``.
pub type Text = Vec<u8>;
/// Type alias for a text slice, i.e. ``&[u8]``.
pub type TextSlice<'a> = &'a [u8];

/// Maximum size of temporary buffer used for reading indexed FASTA files.
const MAX_FASTA_BUFFER_SIZE: usize = 512;

/// Trait for FASTA readers.
pub trait FastaRead {
    fn read(&mut self, record: &mut Record) -> io::Result<()>;
}

/// A FASTA reader.
#[derive(Debug)]
pub struct Reader<B> {
    reader: B,
    line: String,
}

impl Reader<io::BufReader<fs::File>> {
    /// Read FASTA from given file path.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> anyhow::Result<Self> {
        fs::File::open(&path)
            .map(Reader::new)
            .with_context(|| format!("Failed to read fasta from {:#?}", path))
    }

    /// Read FASTA from give file path and a capacity
    pub fn from_file_with_capacity<P: AsRef<Path> + std::fmt::Debug>(
        capacity: usize,
        path: P,
    ) -> anyhow::Result<Self> {
        fs::File::open(&path)
            .map(|file| Reader::with_capacity(capacity, file))
            .with_context(|| format!("Failed to read fasta from {:#?}", path))
    }
}

impl<R> Reader<io::BufReader<R>>
where
    R: io::Read,
{
    /// Create a new Fasta reader given an instance of `io::Read`.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # use bio::io::fasta::Reader;
    /// # fn main() {
    /// # const fasta_file: &'static [u8] = b">id desc
    /// # AAAA
    /// # ";
    /// let reader = Reader::new(fasta_file);
    /// # }
    /// ```
    pub fn new(reader: R) -> Self {
        Reader {
            reader: io::BufReader::new(reader),
            line: String::new(),
        }
    }

    /// Create a new Fasta reader given a capacity and an instance of `io::Read`.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # use bio::io::fasta::Reader;
    /// # fn main() {
    /// # const fasta_file: &'static [u8] = b">id desc
    /// # AAAA
    /// # ";
    /// let reader = Reader::with_capacity(16384, fasta_file);
    /// # }
    /// ```
    pub fn with_capacity(capacity: usize, reader: R) -> Self {
        Reader {
            reader: io::BufReader::with_capacity(capacity, reader),
            line: String::new(),
        }
    }
}

impl<B> Reader<B>
where
    B: io::BufRead,
{
    /// Create a new Fasta reader with an object that implements `io::BufRead`.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # use bio::io::fasta::Reader;
    /// # fn main() {
    /// # const fasta_file: &'static [u8] = b">id desc
    /// # AAAA
    /// # ";
    /// let buffer = io::BufReader::with_capacity(16384, fasta_file);
    /// let reader = Reader::from_bufread(buffer);
    /// # }
    /// ```
    pub fn from_bufread(bufreader: B) -> Self {
        Reader {
            reader: bufreader,
            line: String::new(),
        }
    }

    /// Return an iterator over the records of this Fasta file.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # use bio::io::fasta::Reader;
    /// # use bio::io::fasta::Record;
    /// # fn main() {
    /// # const fasta_file: &'static [u8] = b">id desc
    /// # AAAA
    /// # ";
    /// # let reader = Reader::new(fasta_file);
    /// for record in reader.records() {
    ///     let record = record.unwrap();
    ///     assert_eq!(record.id(), "id");
    ///     assert_eq!(record.desc().unwrap(), "desc");
    ///     assert_eq!(record.seq().to_vec(), b"AAAA");
    /// }
    /// # }
    /// ```
    pub fn records(self) -> Records<B> {
        Records {
            reader: self,
            error_has_occured: false,
        }
    }
}

impl<B> FastaRead for Reader<B>
where
    B: io::BufRead,
{
    /// Read the next FASTA record into the given `Record`.
    /// An empty record indicates that no more records can be read.
    ///
    /// Use this method when you want to read records as fast as
    /// possible because it allows the reuse of a `Record` allocation.
    ///
    /// The [records](Reader::records) iterator provides a more ergonomic
    /// approach to accessing FASTA records.
    ///
    /// # Errors
    ///
    /// This function will return an error if the record is incomplete,
    /// syntax is violated or any form of I/O error is encountered.
    ///
    /// # Example
    ///
    /// ```rust
    /// use bio::io::fasta::Record;
    /// use bio::io::fasta::{FastaRead, Reader};
    ///
    /// const fasta_file: &'static [u8] = b">id desc
    /// AAAA
    /// ";
    /// let mut reader = Reader::new(fasta_file);
    /// let mut record = Record::new();
    ///
    /// // Check for errors parsing the record
    /// reader
    ///     .read(&mut record)
    ///     .expect("fasta reader: got an io::Error or could not read_line()");
    ///
    /// assert_eq!(record.id(), "id");
    /// assert_eq!(record.desc().unwrap(), "desc");
    /// assert_eq!(record.seq().to_vec(), b"AAAA");
    /// ```
    fn read(&mut self, record: &mut Record) -> io::Result<()> {
        record.clear();
        if self.line.is_empty() {
            self.reader.read_line(&mut self.line)?;
            if self.line.is_empty() {
                return Ok(());
            }
        }

        if !self.line.starts_with('>') {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Expected > at record start.",
            ));
        }
        let mut header_fields = self.line[1..].trim_end().splitn(2, char::is_whitespace);
        record.id = header_fields.next().map(|s| s.to_owned()).unwrap();
        record.desc = header_fields.next().map(|s| s.to_owned());
        loop {
            self.line.clear();
            self.reader.read_line(&mut self.line)?;
            if self.line.is_empty() || self.line.starts_with('>') {
                break;
            }
            record.seq.push_str(self.line.trim_end());
        }

        Ok(())
    }
}


/// A Fasta writer.
#[derive(Debug)]
pub struct Writer<W: io::Write> {
    writer: io::BufWriter<W>,
}

impl Writer<fs::File> {
    /// Write to the given file path.
    #[allow(clippy::wrong_self_convention)]
    pub fn to_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::create(path).map(Writer::new)
    }

    /// Write to the given file path and a buffer capacity
    pub fn to_file_with_capacity<P: AsRef<Path>>(capacity: usize, path: P) -> io::Result<Self> {
        fs::File::create(path).map(|file| Writer::with_capacity(capacity, file))
    }
}

impl<W: io::Write> Writer<W> {
    /// Create a new Fasta writer.
    pub fn new(writer: W) -> Self {
        Writer {
            writer: io::BufWriter::new(writer),
        }
    }

    /// Create a new Fasta writer with a capacity of write buffer
    pub fn with_capacity(capacity: usize, writer: W) -> Self {
        Writer {
            writer: io::BufWriter::with_capacity(capacity, writer),
        }
    }

    /// Create a new Fasta writer with a given BufWriter
    pub fn from_bufwriter(bufwriter: io::BufWriter<W>) -> Self {
        Writer { writer: bufwriter }
    }

    /// Directly write a [`fasta::Record`](struct.Record.html).
    ///
    /// # Errors
    /// If there is an issue writing to the `Writer`.
    ///
    /// # Examples
    /// ```rust
    /// use bio::io::fasta::{Record, Writer};
    /// use std::fs;
    /// use std::io;
    /// use std::path::Path;
    ///
    /// let path = Path::new("test.fa");
    /// let file = fs::File::create(path).unwrap();
    /// {
    ///     let handle = io::BufWriter::new(file);
    ///     let mut writer = Writer::new(handle);
    ///     let record = Record::with_attrs("id", Some("desc"), b"ACGT");
    ///
    ///     let write_result = writer.write_record(&record);
    ///     assert!(write_result.is_ok());
    /// }
    ///
    /// let actual = fs::read_to_string(path).unwrap();
    /// let expected = ">id desc\nACGT\n";
    ///
    /// assert!(fs::remove_file(path).is_ok());
    /// assert_eq!(actual, expected)
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        self.write(record.id(), record.desc(), record.seq())
    }

    /// Write a Fasta record with given id, optional description and sequence.
    pub fn write(&mut self, id: &str, desc: Option<&str>, seq: TextSlice<'_>) -> io::Result<()> {
        self.writer.write_all(b">")?;
        self.writer.write_all(id.as_bytes())?;
        if let Some(desc) = desc {
            self.writer.write_all(b" ")?;
            self.writer.write_all(desc.as_bytes())?;
        }
        self.writer.write_all(b"\n")?;
        self.writer.write_all(seq)?;
        self.writer.write_all(b"\n")?;

        Ok(())
    }

    /// Flush the writer, ensuring that everything is written.
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

/// A FASTA record.
#[derive(Default, Clone, Debug, Serialize, Deserialize)]
pub struct Record {
    id: String,
    desc: Option<String>,
    seq: String,
}

impl Record {
    /// Create a new instance.
    pub fn new() -> Self {
        Record {
            id: String::new(),
            desc: None,
            seq: String::new(),
        }
    }

    /// Create a new `Record` from given attributes.
    ///
    /// # Examples
    /// ```rust
    /// use bio::io::fasta::Record;
    ///
    /// let read_id = "read1";
    /// let description = Some("sampleid=foobar");
    /// let sequence = b"ACGT";
    /// let record = Record::with_attrs(read_id, description, sequence);
    ///
    /// assert_eq!(">read1 sampleid=foobar\nACGT\n", record.to_string())
    /// ```
    pub fn with_attrs(id: &str, desc: Option<&str>, seq: TextSlice<'_>) -> Self {
        let desc = desc.map(|desc| desc.to_owned());
        Record {
            id: id.to_owned(),
            desc,
            seq: String::from_utf8(seq.to_vec()).unwrap(),
        }
    }

    /// Check if record is empty.
    pub fn is_empty(&self) -> bool {
        self.id.is_empty() && self.desc.is_none() && self.seq.is_empty()
    }

    /// Check validity of Fasta record.
    pub fn check(&self) -> Result<(), &str> {
        if self.id().is_empty() {
            return Err("Expecting id for Fasta record.");
        }
        if !self.seq.is_ascii() {
            return Err("Non-ascii character found in sequence.");
        }

        Ok(())
    }

    /// Return the id of the record.
    pub fn id(&self) -> &str {
        self.id.as_ref()
    }

    /// Return descriptions if present.
    pub fn desc(&self) -> Option<&str> {
        match self.desc.as_ref() {
            Some(desc) => Some(desc),
            None => None,
        }
    }

    /// Return the sequence of the record.
    pub fn seq(&self) -> TextSlice<'_> {
        self.seq.as_bytes()
    }

    /// Clear the record.
    fn clear(&mut self) {
        self.id.clear();
        self.desc = None;
        self.seq.clear();
    }
}

impl fmt::Display for Record {
    /// Allows for using `Record` in a given formatter `f`. In general this is for
    /// creating a `String` representation of a `Record` and, optionally, writing it to
    /// a file.
    ///
    /// # Errors
    /// Returns [`std::fmt::Error`](https://doc.rust-lang.org/std/fmt/struct.Error.html)
    /// if there is an issue formatting to the stream.
    ///
    /// # Examples
    ///
    /// Read in a Fasta `Record` and create a `String` representation of it.
    ///
    /// ```rust
    /// use bio::io::fasta::Reader;
    /// use std::fmt::Write;
    /// // create a "fake" fasta file
    /// let fasta: &'static [u8] = b">id comment1 comment2\nACGT\n";
    /// let mut records = Reader::new(fasta).records().map(|r| r.unwrap());
    /// let record = records.next().unwrap();
    ///
    /// let mut actual = String::new();
    /// // populate `actual` with a string representation of our record
    /// write!(actual, "{}", record).unwrap();
    ///
    /// let expected = std::str::from_utf8(fasta).unwrap();
    ///
    /// assert_eq!(actual, expected)
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let header = match self.desc() {
            Some(d) => format!("{} {}", self.id().to_owned(), d),
            None => self.id().to_owned(),
        };
        write!(
            f,
            ">{}\n{}\n",
            header,
            std::str::from_utf8(self.seq()).unwrap(),
        )
    }
}

/// An iterator over the records of a Fasta file.
pub struct Records<B>
where
    B: io::BufRead,
{
    reader: Reader<B>,
    error_has_occured: bool,
}

impl<B> Iterator for Records<B>
where
    B: io::BufRead,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<io::Result<Record>> {
        if self.error_has_occured {
            None
        } else {
            let mut record = Record::new();
            match self.reader.read(&mut record) {
                Ok(()) if record.is_empty() => None,
                Ok(()) => Some(Ok(record)),
                Err(err) => {
                    self.error_has_occured = true;
                    Some(Err(err))
                }
            }
        }
    }
}
