import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_protein_sequence(sequence):
    """
    Identifies a protein sequence using NCBI BLAST.
    
    Args:
        sequence (str): The amino acid sequence to identify.
    """
    try:
        # Clean up the sequence by removing spaces and newlines
        query_sequence = "".join(sequence.split())
        
        print("Performing BLAST search against Swiss-Prot database...")
        print("This may take a moment.")
        
        # Perform the BLAST search online
        # Using "blastp" for protein-protein comparison
        # Using "swissprot" for a curated database
        result_handle = NCBIWWW.qblast("blastp", "swissprot", query_sequence)
        
        # Parse the BLAST results
        blast_record = NCBIXML.read(result_handle)
        
        # Check if any alignments were found
        if not blast_record.alignments:
            print("No significant matches found for the provided sequence.")
            return

        # Get the top hit (the first alignment)
        top_alignment = blast_record.alignments[0]
        
        # The title contains the protein description. 
        # e.g., >sp|P27348.3|IDH1_HUMAN Isocitrate dehydrogenase [NADP] cytoplasmic OS=...
        # We parse this title to get the common protein name.
        protein_title = top_alignment.title
        protein_description = protein_title.split(' OS=')[0].split(' ', 1)[1]

        print("\n--- BLAST Search Result ---")
        print(f"Top Match Title: {protein_title}")
        print(f"Identified Protein: {protein_description}")
        print(f"E-value (Expect value): {top_alignment.hsps[0].expect}")
        print(f"Identity: {top_alignment.hsps[0].identities}/{top_alignment.hsps[0].align_length} ({100*top_alignment.hsps[0].identities/top_alignment.hsps[0].align_length:.2f}%)")
        
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        print("Please check your internet connection and ensure 'biopython' is installed ('pip install biopython').", file=sys.stderr)

if __name__ == "__main__":
    # The amino acid sequence provided by the user
    input_sequence = """
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
    identify_protein_sequence(input_sequence)