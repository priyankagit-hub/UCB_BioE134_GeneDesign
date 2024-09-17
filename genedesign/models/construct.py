from dataclasses import dataclass
from typing import List
from .transcript import Transcript  # Assuming Transcript is defined in transcript.py

@dataclass(frozen=True)
class Construct:
    """
    Encodes a genetic construct described in terms of a single operon
    encoding multiple cocistronic transcripts.
    """
    mRNAs: List[Transcript]
    promoter: str
    terminator: str
