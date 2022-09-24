# Base error classes
class __BaseFormatError(Exception):
    """Captures any format related issues"""

class __BaseIterationError(Exception):
    pass

# Specific Exceptions
class FastqFormatError(__BaseFormatError):
    """Raised when any format related issues with Fastq files"""

    
class SearchExceededError(__BaseIterationError):
    """Raised if any searching algorithm exceeds max search iteration"""
