import sys
import io

# The biopython library is required for this script.
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Error: The 'biopython' library is required. Please install it using 'pip install biopython' and try again.", file=sys.stderr)
    sys.exit(1)

def identify_protein_from_sequence(sequence_data):
    """
    Identifies a protein by performing a BLAST search with its amino acid sequence.

    Args:
        sequence_data (str): A string containing the amino acid sequence, possibly with whitespace and newlines.

    Returns:
        str: The name of the most likely protein, or an error message.
    """
    # Clean up the sequence by removing whitespace and newlines
    protein_sequence = "".join(sequence_data.split())

    print("Performing BLAST search against the NCBI database...")
    print("This may take a moment, please be patient.")

    result_handle = None
    try:
        # Use qblast to send the sequence to the NCBI BLAST server
        # program='blastp' for protein-protein blast
        # database='nr' for the non-redundant protein database
        # Using StringIO to handle the result in memory
        result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence)
        
        # We use read() as we expect only one query result for a single sequence
        blast_record = NCBIXML.read(result_handle)

        # Check if any alignments were found
        if blast_record.alignments:
            # The first alignment object contains the best hit
            top_alignment = blast_record.alignments[0]
            # The title of the hit sequence contains the protein name and organism
            protein_name = top_alignment.title
            
            # The HSP (High-scoring Segment Pair) contains statistics
            top_hsp = top_alignment.hsps[0]
            e_value = top_hsp.expect
            identity_percent = (top_hsp.identities / top_hsp.align_length) * 100

            print("\n--- Search Complete ---")
            print(f"Top Match Found: {protein_name}")
            print(f"Significance (E-value): {e_value:.2e}")
            print(f"Identity: {identity_percent:.2f}%")
            
            # Extract the core protein name for the final answer
            # Example title: "gi|...|ref|XP_...| argininosuccinate synthase 1 [Homo sapiens]"
            final_answer = protein_name.split(']')[0].split('|')[-1].strip()
            if '[' in final_answer: # handle cases where the name is before the bracket
                final_answer = final_answer.split('[')[0].strip()
            return final_answer

        else:
            return "No significant matches found for the given sequence."

    except Exception as e:
        return f"An error occurred during the BLAST search: {e}\nPlease check your internet connection."
    finally:
        if result_handle:
            result_handle.close()

if __name__ == '__main__':
    # The amino acid sequence provided by the user
    sequence_input = """
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
    
    result = identify_protein_from_sequence(sequence_input)
    print(f"\nConclusion: The sequence represents a form of '{result}'.")
