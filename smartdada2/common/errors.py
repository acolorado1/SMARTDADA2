# Base error classes
class __BaseFormatError(Exception):
    """Captures any format related issues"""


# Specific Exceptions
class FastqFormatError(__BaseFormatError):
    """Raised when any format related issues with Fastq files"""
