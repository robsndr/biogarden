
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
    pub fn new(reader: R) -> Self {
        Reader {
            reader: io::BufReader::new(reader),
            line: String::new(),
        }
    }

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

    pub fn from_bufread(bufreader: B) -> Self {
        Reader {
            reader: bufreader,
            line: String::new(),
        }
    }

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
