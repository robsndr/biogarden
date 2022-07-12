use crate::ds::graph::GraphErr;
use std::error::Error;
use std::fmt;


#[derive(Debug)] 
/// Error variants
pub enum BioError {
    InvalidInputSize,
    InvalidArgumentRange,
    ItemNotFound,
    TypeConversionError, 
    GraphError(GraphErr),
}

impl fmt::Display for BioError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            BioError::InvalidInputSize => write!(f, "Provided inputs have invalid size!"),
            BioError::InvalidArgumentRange => write!(f, "The provided has is within an unsupported range!"),
            BioError::ItemNotFound => write!(f, "The requested item does not exist!"),
            BioError::TypeConversionError => write!(f, "The requested type conversion resulted in an error!"),
            BioError::GraphError(ref source) => write!(f, "An error occurred during graph processing! {}", source),
        }
    }
}

impl Error for BioError {
    fn source(&self) -> Option<&(dyn Error + 'static)>{
        match self {
            BioError::InvalidInputSize => None,
            BioError::InvalidArgumentRange => None, 
            BioError::ItemNotFound => None,
            BioError::TypeConversionError => None,
            BioError::GraphError(ref source) => Some(source)
        }
    }
}

impl From<GraphErr> for BioError {
    fn from(cause: GraphErr) -> BioError {
        BioError::GraphError(cause)
    }
}

pub type Result<T> = std::result::Result<T, BioError>;
