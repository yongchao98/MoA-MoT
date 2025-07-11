# First, ensure you have BioPython installed:
# pip install biopython

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_protein(sequence):
    """
    Identifies a protein sequence using NCBI BLAST.
    """
    # Remove any whitespace or newlines from the sequence
    protein_sequence = "".join(sequence.split())
    
    print("Performing BLAST search... This may take a few minutes.")
    
    try:
        # Use qblast to perform the search against the non-redundant (nr) protein database
        result_handle = NCBIWWW.qblast(
            program="blastp", 
            database="nr", 
            sequence=protein_sequence
        )
        
        # Parse the XML result
        blast_record = NCBIXML.read(result_handle)
        
        print("\n--- BLAST Search Results ---")
        if not blast_record.alignments:
            print("No significant matches found.")
            return

        # Print the top 5 results
        limit = 5
        print(f"Showing top {min(len(blast_record.alignments), limit)} matches:\n")
        
        for i, alignment in enumerate(blast_record.alignments):
            if i >= limit:
                break
            for hsp in alignment.hsps:
                print(f"**** Match {i+1} ****")
                print(f"Description: {alignment.title}")
                print(f"Length: {alignment.length}")
                print(f"E-value: {hsp.expect}")
                print(f"Score: {hsp.score}")
                print("-" * 20)
                # We only need to show the first high-scoring pair (HSP) for each alignment
                break
                
    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please check your internet connection and the sequence format.")
    finally:
        if 'result_handle' in locals():
            result_handle.close()

# The amino acid sequence provided by the user
sequence = """
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

identify_protein(sequence)