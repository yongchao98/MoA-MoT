# First, ensure you have biopython installed:
# pip install biopython

from Bio.Blast import NCBIWWW, NCBIXML

def identify_protein(sequence_data):
    """
    Identifies a protein by performing a BLAST search with its amino acid sequence.
    """
    # Clean up the sequence by removing spaces and newlines
    cleaned_sequence = "".join(sequence_data.split())

    print(f"Identifying protein for the sequence: {cleaned_sequence[:30]}...")
    print("Performing BLAST search against the NCBI 'nr' database. This may take a moment...")

    try:
        # Perform the BLAST search using NCBI's web service
        # 'blastp' is for protein-protein comparison
        # 'nr' is the non-redundant protein database
        result_handle = NCBIWWW.qblast("blastp", "nr", cleaned_sequence)

        # Parse the XML results
        blast_record = NCBIXML.read(result_handle)

        # Check if any results were found
        if blast_record.alignments:
            # The first alignment is the best hit
            top_alignment = blast_record.alignments[0]
            protein_title = top_alignment.title
            top_hsp = top_alignment.hsps[0] # High-scoring Segment Pair

            print("\n--- Search Complete ---")
            print("The most likely protein identity is:")
            # The title often contains the protein name and organism in brackets
            # e.g., ">sp|P11712.3| CP2C9_HUMAN RecName: Full=Cytochrome P450 2C9... [Homo sapiens]"
            print(protein_title)
            print(f"\nE-value (Expect value): {top_hsp.expect}")
            print("The E-value represents the number of hits one can 'expect' to see by chance. A lower value indicates a more significant match.")

        else:
            print("No significant matches found for the provided sequence.")

    except Exception as e:
        print(f"An error occurred during the BLAST search: {e}")
        print("This could be due to a temporary issue with the NCBI servers or your internet connection.")
    finally:
        # Clean up the result handle
        if 'result_handle' in locals() and result_handle:
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

if __name__ == '__main__':
    identify_protein(sequence)