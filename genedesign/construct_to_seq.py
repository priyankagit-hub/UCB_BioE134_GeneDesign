from models.construct import Construct
from transcript_to_seq import transcript_to_seq

def construct_to_seq(construct: Construct) -> str:
    """
    Converts a Construct object into its full DNA sequence by concatenating the promoter,
    the sequences of the mRNAs, and the terminator.
    
    Parameters:
        construct (Construct): The construct object containing mRNAs, promoter, and terminator.
    
    Returns:
        str: The full DNA sequence of the construct.
    """
    # Start with the promoter
    out = [construct.promoter]
    
    # Append the sequence of each mRNA by calling transcript_to_seq on each one
    out.extend(transcript_to_seq(mrna) for mrna in construct.mRNAs)
    
    # Finally, append the terminator
    out.append(construct.terminator)
    
    # Join the list into a single string and return
    return ''.join(out)
