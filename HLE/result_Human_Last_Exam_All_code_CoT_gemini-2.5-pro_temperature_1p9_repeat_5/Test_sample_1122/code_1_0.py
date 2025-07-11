# This script requires the 'biopython' library.
# If you don't have it installed, run this command in your terminal or command prompt:
# pip install biopython

import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_protein_from_sequence(sequence):
    """
    Identifies a protein by performing a BLAST search against the NCBI 'nr' database.
    It prints the description of the top hit.
    """
    try:
        print("Submitting BLAST search to NCBI. This may take a minute...")
        
        # Use NCBI's online BLAST service (qblast)
        # 'blastp' is for protein-protein search
        # 'nr' is the non-redundant protein database
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
        
        # The result from NCBI is in XML format. We need to parse it.
        # We use read() because we submitted a single query and expect a single result.
        blast_record = NCBIXML.read(result_handle)
        
        # Check if the search yielded any results (alignments)
        if blast_record.alignments:
            # The first alignment in the record is the best hit
            first_alignment = blast_record.alignments[0]
            
            print("\n--- Top BLAST Hit ---")
            # The title of the alignment contains the protein's description
            print(f"Description: {first_alignment.title}")
            print(f"Length: {first_alignment.length}")
            # The first High-scoring Segment Pair (HSP) contains alignment details
            first_hsp = first_alignment.hsps[0]
            print(f"E-value (Expect): {first_hsp.expect}")
            print(f"Score: {first_hsp.score}")
            print("----------------------")
        else:
            print("No significant matches were found for the given sequence.")
            
    except ImportError:
        print("ERROR: The 'biopython' library is required but not installed.")
        print("Please install it by running: pip install biopython")
    except Exception as e:
        print(f"An error occurred during the BLAST search: {e}")
        print("This could be due to a network issue or a problem with the NCBI servers.")
    finally:
        # It's good practice to close the handle
        if 'result_handle' in locals() and result_handle is not None:
            result_handle.close()

if __name__ == "__main__":
    # The amino acid sequence provided by the user
    amino_acid_sequence = (
        "MAEQVALSRT QVCGILREEL FQGDAFHQSD THIFIIMGAS GDLAKKKIYP "
        "TIWWLFRDGL LPENTFIVGY ARSRLTVADI RKQSEPFFKA TPEEKLKLED "
        "FFARNSYVAG QYDDAASYQR LNSHMNALHL GSQANRLFYL ALPPTVYEAV "
        "TKNIHESCMQ IGWNRIIVEK PFGRDLQSSD RLSNHISSLF REDQIYRIDH "
        "YLGKEMVQNL MVLRFANRIF GPIWNRDNIA CVILTFKEPF GTEGRGGYFD "
        "EFGIIRDVMQ NHLLQMLCLV AMEKPASTNS DDVRDEKVKV LKCISEVQAN "
        "NVVLGQYVGN PDGEGEATKG YLDDPTVPRG STTATFAAVV LYVENERWDG "
        "VPFILRCGKA LNERKAEVRL QFHDVAGDIF HQQCKRNELV IRVQPNEAVY "
        "TKMMTKKPGM FFNPEESELD LTYGNRYKNV KLPDAYERLI LDVFCGSQMH "
        "FVRSDELREA WRIFTPLLHQ IELEKPKPIP YIYGSRGPTE ADELMKRVGF "
        "QYEGTYKWVN PHKL"
    )
    
    # Clean the sequence by removing spaces and newline characters
    cleaned_sequence = "".join(amino_acid_sequence.split())
    
    identify_protein_from_sequence(cleaned_sequence)
