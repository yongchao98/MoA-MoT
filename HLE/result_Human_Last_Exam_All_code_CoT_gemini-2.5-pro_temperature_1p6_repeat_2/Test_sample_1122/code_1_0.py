import io
from xml.etree import ElementTree

try:
    # Biopython is required for this script.
    from Bio.Blast import NCBIWWW
except ImportError:
    print("Biopython library not found.")
    print("Please install it using: pip install biopython")
    exit()

def identify_protein(sequence: str):
    """
    Identifies a protein by performing a BLAST search against the NCBI database.
    
    Args:
        sequence (str): The amino acid sequence of the protein.
    """
    # Remove any whitespace or newlines from the sequence string
    formatted_sequence = "".join(sequence.split())
    
    print("Performing BLAST search with the provided sequence...")
    print("This may take a moment, depending on the NCBI server load.")
    
    try:
        # Use NCBI's qblast function. 
        # program="blastp" for protein-protein comparison
        # database="nr" for the non-redundant protein database
        # sequence=the amino acid sequence
        result_handle = NCBIWWW.qblast("blastp", "nr", formatted_sequence)

        # The result from NCBI is in XML format. We read it as a string.
        blast_xml_result = result_handle.read()
        result_handle.close()

        # To avoid printing the entire raw XML, we parse it to find the top hit.
        # We can treat the XML string as a file using io.StringIO
        xml_io = io.StringIO(blast_xml_result)
        tree = ElementTree.parse(xml_io)
        root = tree.getroot()
        
        # The path in the XML to the description of the first hit is:
        # BlastOutput_iterations/Iteration/Iteration_hits/Hit/Hit_def
        first_hit_description = root.find('.//Hit_def')

        if first_hit_description is not None and first_hit_description.text:
            print("\n--- BLAST Search Result ---")
            print("The most likely protein is:")
            print(first_hit_description.text)
            # The format is requested for the final answer block
            # but we can return it from the function too.
            return first_hit_description.text
        else:
            print("\nCould not find a significant match for the sequence.")
            # Check for potential errors in the XML response
            if root.find('.//BlastOutput_error'):
                error_msg = root.find('.//BlastOutput_error').text
                print(f"An error was reported by the server: {error_msg}")
            return "No significant match found."

    except Exception as e:
        print(f"\nAn error occurred during the BLAST search: {e}")
        print("This may be due to a network issue or a problem with the NCBI BLAST server.")
        return "Error during search"

if __name__ == "__main__":
    # The amino acid sequence from the user
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
    
    protein_identity = identify_protein(amino_acid_sequence)
