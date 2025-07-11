# First, ensure you have the biopython library installed:
# pip install biopython

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import io

def identify_protein_from_sequence(sequence):
    """
    Identifies a protein by performing a BLAST search with its amino acid sequence.
    """
    # Step 1: Clean up the sequence string
    cleaned_sequence = "".join(sequence.strip().split())
    
    print("Performing BLAST search for the protein sequence.")
    print("This may take a few moments...")
    
    try:
        # Step 2: Perform the BLAST search using NCBI's online service
        # Program: 'blastp' (protein-protein)
        # Database: 'nr' (non-redundant protein sequences)
        result_handle = NCBIWWW.qblast("blastp", "nr", cleaned_sequence)
        
        # Step 3: Parse the XML results
        # We read the result into a string and then parse it to avoid handle issues
        blast_xml_string = result_handle.read()
        result_handle.close()
        
        # Use io.StringIO to treat the XML string as a file for parsing
        blast_record = NCBIXML.read(io.StringIO(blast_xml_string))
        
        # Step 4: Display the top result
        if blast_record.descriptions:
            print("\n--- Top Match Found ---")
            top_description = blast_record.descriptions[0]
            print(f"Description: {top_description.title}")
            print(f"Score: {top_description.score}")
            print(f"E-value: {top_description.e}")
            # The full name is often available as 'Full=...' in the description
            try:
                full_name_start = top_description.title.index('Full=') + 5
                full_name_end = top_description.title.index(';', full_name_start)
                protein_name = top_description.title[full_name_start:full_name_end]
            except ValueError:
                # Fallback if 'Full=' is not found
                protein_name = top_description.title.split(' [')[0]

            print(f"\nThis sequence most likely represents the protein: {protein_name}")
            return protein_name
        else:
            print("\nNo significant matches were found for the provided sequence.")
            return None

    except Exception as e:
        print(f"\nAn error occurred: {e}")
        print("Please check your internet connection or try again later.")
        return None

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

# Run the identification function
identify_protein_from_sequence(input_sequence)