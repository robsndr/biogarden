use std::io::Error as IoError;
use std::str::Utf8Error;
use std::convert::From;
use std::error::Error;
use std::fmt;

#[derive(Debug)] // Allow the use of "{:?}" format specifier
pub enum BioError {
    InvalidInputSize,
}

// Allow the use of "{}" format specifier
impl fmt::Display for BioError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            BioError::InvalidInputSize => write!(f, "Provided inputs have invalid size!"),
        }
    }
}

// Allow this type to be treated like an error
impl Error for BioError {
    fn source(&self) -> Option<&(dyn Error + 'static)>{
        match *self {
            BioError::InvalidInputSize => None,
        }
    }
}

// // Support converting system errors into our custom error.
// // This trait is used in `try!`.
// impl From<IoError> for CustomError {
//     fn from(cause: IoError) -> CustomError {
//         CustomError::Io(cause)
//     }
// }
// impl From<Utf8Error> for CustomError {
//     fn from(cause: Utf8Error) -> CustomError {
//         CustomError::Utf8(cause)
//     }
// }

pub type Result<T> = std::result::Result<T, BioError>;
