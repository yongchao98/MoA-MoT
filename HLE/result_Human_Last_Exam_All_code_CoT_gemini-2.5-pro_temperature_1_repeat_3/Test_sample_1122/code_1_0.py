# First, you need to install the BioPython library if you haven't already:
# pip install biopython

import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_protein_from_sequence(sequence):
    """
    Identifies a protein from its amino acid sequence using NCBI BLAST.
    Requires an active internet connection.
    """
    # Clean up the sequence by removing spaces and newlines
    cleaned_sequence = "".join(sequence.split())
    
    print("Performing BLAST search against the NCBI database...")
    print("This may take a moment, please be patient.")
    
    try:
        # Use qblast to perform the search online.
        # program: blastp (protein-protein)
        # database: nr (non-redundant protein database)
        # sequence: our cleaned sequence
        result_handle = NCBIWWW.qblast("blastp", "nr", cleaned_sequence)

        # Parse the XML results
        # NCBIXML.read() is used for a single query result
        blast_record = NCBIXML.read(result_handle)

        # Check if any results were found
        if blast_record.alignments:
            # The first alignment is the best hit
            top_hit = blast_record.alignments[0]
            print("\n--- Top BLAST Hit Found ---")
            print(f"Title: {top_hit.title}")
            print(f"Length: {top_hit.length}")
            # The hsps object contains alignment details
            # We'll print the E-value (Expect value) of the top hit
            # A lower E-value indicates a more significant match.
            top_hsp = top_hit.hsps[0]
            print(f"E-value: {top_hsp.expect}")
        else:
            print("\nNo significant match was found for the provided sequence.")

    except Exception as e:
        print(f"\nAn error occurred during the BLAST search: {e}", file=sys.stderr)
        print("Please check your internet connection and ensure the sequence is a valid amino acid sequence.", file=sys.stderr)
    finally:
        # It's good practice to close the handle
        if 'result_handle' in locals() and result_handle:
            result_handle.close()


# The amino acid sequence provided by the user
amino_acid_sequence = """
MAEQVALSRT QVCGILREEL FQGDAFHQSD THIFIIMGAS GDLAKKKIYP 
TIWWLFRDGL LPENTFIVGY ARSRLTVADI RKQSEPFFKA TPEEKLKLED 
FFARNSYVAG QYDDAASYQR LNSHMNALHL GSQANRLFYL ALPPTVYEAV 
TKNIHESCMQ IGWNRIIVEK PFGRDLQSSD RLSNHISSLF REDQIYRIDH 
YLGKEMVQNL MVLRFANRIF GPIWNRDNIA CVILTFKEPF GTEGRGGYFD 
EFGIIRDVMQ NHLLQMLCLV AMEKPASTNS DDVRDEKVKV LKCISEVQAN 
NVVLGQYVGN PDGEGEATKG YLDDPTVPRG STTATFAAVV LYVENERWDG 
VPFILRCGKA LNERKAEVRL QFHDVAGDIF HQQCKRNELV IRVQPNEAVY 
TKMMTKKPGM FFNPEESELD LTYGNRYKNV KLPDAYERLI LDVFCGSQMH 
FVRSDELREA WRIFTPLLHQ IELEKPKPIP YIYGSRGPTE ADELMKRVGF 
QYEGTYKWVN PHKL
"""

# Run the function with the user's sequence
identify_protein_from_sequence(amino_acid_sequence)