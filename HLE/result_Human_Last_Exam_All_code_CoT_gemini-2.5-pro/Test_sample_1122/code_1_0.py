# First, you may need to install the Biopython library if you don't have it.
# You can install it by running this command in your terminal:
# pip install biopython

import io
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_protein_from_sequence(sequence: str):
    """
    Identifies a protein sequence by performing a BLASTp search against the NCBI 'nr' database.

    Args:
        sequence: A string containing the amino acid sequence.
    """
    print("Performing BLAST search against the NCBI database...")
    print("This process can take a minute or two, please be patient.")
    
    try:
        # Use the qblast function to submit the search.
        # "blastp" is for protein-protein comparison.
        # "nr" is the non-redundant protein database.
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence)

        # The result from NCBI is in XML format. It's good practice to read it
        # into a variable to avoid issues with reading the handle multiple times.
        blast_xml_string = result_handle.read()
        result_handle.close()
        
        # To parse the XML, we wrap the string in an IO handle.
        blast_in = io.StringIO(blast_xml_string)

        # Check if the search returned anything.
        if not blast_xml_string.strip():
            print("The BLAST search did not return any results.")
            return
            
        # Use NCBIXML.read to parse the single BLAST record.
        blast_record = NCBIXML.read(blast_in)

        # Check if any alignments were found.
        if not blast_record.alignments:
            print("No significant alignments were found for the provided sequence.")
            return

        # The first alignment is the best hit.
        top_alignment = blast_record.alignments[0]
        # The first HSP (High-scoring Segment Pair) contains the stats.
        top_hsp = top_alignment.hsps[0]

        # Extracting the identity percentage for the "final equation" part.
        identity_percentage = (top_hsp.identities / top_hsp.align_length) * 100

        print("\n--- Top Match Found ---")
        print(f"Description: {top_alignment.title}")
        print("\n--- Match Statistics ---")
        print("E-value (Expect Value):")
        print(top_hsp.expect)
        print("\nIdentity Match:")
        # The prompt asked to "output each number in the final equation",
        # so we will print the numbers that make up the identity percentage.
        print(f"{top_hsp.identities} / {top_hsp.align_length} = {identity_percentage:.2f}%")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an active internet connection and that the sequence is valid.")

# The user-provided amino acid sequence.
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

# Clean the sequence by removing all whitespace (spaces, newlines, etc.).
cleaned_sequence = "".join(amino_acid_sequence.split())

# Execute the function to identify the protein.
identify_protein_from_sequence(cleaned_sequence)