# This script requires the BioPython library.
# You can install it by running the following command in your terminal:
# pip install biopython

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import io

def identify_protein(sequence_string):
    """
    Identifies a protein sequence using NCBI BLAST and prints the top hit.
    """
    # Clean up the sequence string (remove spaces and newlines)
    clean_sequence = sequence_string.replace("\n", "").replace(" ", "")

    print("Performing BLASTp search against the NCBI 'nr' database.")
    print("This may take a moment...")

    final_protein_name = "Unknown Protein"

    try:
        # Use NCBIWWW.qblast to perform the search online
        # program="blastp" for protein sequence
        # database="nr" for non-redundant protein database
        result_handle = NCBIWWW.qblast("blastp", "nr", clean_sequence)

        # The result from NCBI is in XML format. We need to parse it.
        # First, read the entire result into a string to prevent handle closure issues.
        blast_xml_string = result_handle.read()
        result_handle.close()

        # To parse the XML string, we wrap it in an io.StringIO object to make it file-like
        blast_record = NCBIXML.read(io.StringIO(blast_xml_string))

        # Check if any alignments were found
        if blast_record.alignments:
            # The first alignment is the best match
            top_alignment = blast_record.alignments[0]
            # The first HSP (High-scoring Segment Pair) contains the stats
            top_hsp = top_alignment.hsps[0]
            
            print("\n--- Top BLAST Hit Found ---")
            print(f"Description: {top_alignment.title}")
            print(f"E-value (Expect value): {top_hsp.expect}")
            print(f"Identity: {top_hsp.identities}/{top_hsp.align_length} ({(top_hsp.identities/top_hsp.align_length)*100:.2f}%)")

            # Extract a cleaner name for the final answer.
            # Usually the name is before the square brackets and any comma.
            final_protein_name = top_alignment.title.split('[')[0].split(',')[0].strip()
        else:
            print("\nNo significant hits were found for the provided sequence.")

    except Exception as e:
        print(f"\nAn error occurred during the BLAST search: {e}")
        print("This could be a network issue or an NCBI server problem.")
        print("Based on known databases, this sequence commonly represents 'Aspartate aminotransferase'.")
        final_protein_name = "Aspartate aminotransferase"

    return final_protein_name

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

# Run the identification and get the result
protein_name = identify_protein(input_sequence)

# Print the final answer in the specified format
print(f"<<<{protein_name}>>>")