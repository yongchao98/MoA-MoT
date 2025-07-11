# First, ensure you have biopython installed:
# pip install biopython

from Bio.Blast import NCBIWWW, NCBIXML
import io

def identify_protein(sequence):
    """
    Identifies a protein sequence using NCBI BLAST.
    """
    # Combine the sequence into a single string without spaces or newlines
    protein_sequence = "".join(sequence.split())

    print(f"Querying NCBI BLAST with sequence:\n{protein_sequence}\n")
    print("This may take a moment...")

    try:
        # Perform the BLAST search
        # program="blastp" for protein sequence
        # database="nr" for non-redundant protein database
        result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence)

        # The result from NCBI is in XML format. We'll read it as a string
        # to use it with NCBIXML.read, which expects a file-like object.
        blast_xml_string = result_handle.read()
        result_handle.close()

        # Use io.StringIO to treat the XML string as a file
        blast_record = NCBIXML.read(io.StringIO(blast_xml_string))

        # Check if any alignments were found
        if not blast_record.alignments:
            print("No significant matches found for the provided sequence.")
            return

        # The first alignment is the best hit
        top_alignment = blast_record.alignments[0]
        top_hsp = top_alignment.hsps[0] # High-scoring Segment Pair

        print("\n--- Top BLAST Hit ---")
        # The title contains the protein name and organism
        print(f"Description: {top_alignment.title}")
        print(f"Accession: {top_alignment.accession}")
        print(f"E-value: {top_hsp.expect}")
        print(f"Identity: {top_hsp.identities}/{top_hsp.align_length} ({100 * top_hsp.identities / top_hsp.align_length:.2f}%)")

        # Extract the main protein name from the title for the final answer
        # The format is often ">accession protein_name [organism]"
        # We'll split by the first space after the accession to get the name
        try:
            # The title often looks like:
            # "sp|P32119.1|KAD2_HUMAN RecName: Full=Aryl-hydrocarbon receptor-interacting protein;..."
            # or "NP_003968.1 aryl-hydrocarbon receptor-interacting protein [Homo sapiens]"
            # We will extract the part after the accession code and before the organism bracket.
            full_title = top_alignment.title
            # Find the first space to skip the accession part
            first_space_index = full_title.find(' ')
            # Find the start of the organism name
            organism_bracket_index = full_title.find('[')
            if first_space_index != -1 and organism_bracket_index != -1:
                protein_name = full_title[first_space_index + 1 : organism_bracket_index].strip()
            else: # Fallback if format is unusual
                 protein_name = " ".join(full_title.split(' ')[1:])

            return protein_name

        except Exception:
            # Fallback if parsing fails
            return "Could not parse the exact protein name from the title."


    except Exception as e:
        print(f"An error occurred during the BLAST search: {e}")
        print("Please check your internet connection and ensure the NCBI services are available.")
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

protein_name = identify_protein(input_sequence)

if protein_name:
    print("\nFinal Answer:")
    # Use triple-quotes for the answer format as requested
    print(f'<<<{protein_name}>>>')
