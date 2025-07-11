# First, you need to install the Biopython library if you haven't already.
# You can do this by running the following command in your terminal or command prompt:
# pip install biopython

import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_protein(sequence):
    """
    Identifies a protein sequence using NCBI's BLAST service.
    """
    print("Performing BLAST search. This may take a minute or two...")
    
    # Use the qblast function to submit the search to NCBI.
    # program='blastp' is for protein-protein comparison.
    # database='nr' is the non-redundant protein database.
    try:
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
    except Exception as e:
        print(f"An error occurred while contacting NCBI: {e}")
        print("Please check your internet connection and try again.")
        return

    # Parse the XML result returned by BLAST
    try:
        blast_record = NCBIXML.read(result_handle)
    except Exception as e:
        print(f"Failed to parse BLAST results: {e}")
        return
        
    # Check if any alignments were found
    if not blast_record.alignments:
        print("No matches found for the provided sequence.")
        return

    # The top hit is the first alignment in the list.
    top_hit = blast_record.alignments[0]
    
    print("\n--- Top BLAST Hit ---")
    print(f"Description: {top_hit.title}")
    
    # An alignment's hsp (High-scoring Segment Pair) contains the detailed
    # alignment information like score, e-value, and identities.
    top_hsp = top_hit.hsps[0]
    
    print(f"Accession ID: {top_hit.accession}")
    print(f"E-value (Expect value): {top_hsp.expect}")
    
    identity_percentage = (top_hsp.identities / top_hsp.align_length) * 100
    print(f"Identity: {top_hsp.identities}/{top_hsp.align_length} ({identity_percentage:.2f}%)")
    print("-----------------------\n")
    
    # Extract just the protein name for the final answer
    # The title is often formatted like: 'gi|...|ref|...| Protein Name [Organism]'
    protein_name_full = top_hit.title.split('|')[-1].strip()
    print(f"The sequence most likely represents: {protein_name_full}")


if __name__ == "__main__":
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

    # Clean up the sequence by removing spaces and newlines
    cleaned_sequence = amino_acid_sequence.replace(" ", "").replace("\n", "")

    identify_protein(cleaned_sequence)
