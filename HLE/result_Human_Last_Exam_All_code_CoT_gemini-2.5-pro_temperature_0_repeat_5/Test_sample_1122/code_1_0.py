import io
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_protein_sequence(sequence_data):
    """
    Identifies a protein sequence using NCBI BLAST.

    Args:
        sequence_data: A string containing the amino acid sequence.
    """
    try:
        # Clean up the sequence by removing spaces and newlines
        cleaned_sequence = "".join(sequence_data.split())

        print("Performing BLAST search against the NCBI database...")
        print("This may take a minute or two, please wait.")

        # Perform the BLAST search using the qblast function
        # 'blastp' is for protein-protein search
        # 'nr' is the non-redundant protein database
        result_handle = NCBIWWW.qblast("blastp", "nr", cleaned_sequence)

        # The result is in XML format. We read it as a string.
        blast_results_xml = result_handle.read()
        result_handle.close()

        # Use io.StringIO to treat the XML string as a file for parsing
        blast_record = NCBIXML.read(io.StringIO(blast_results_xml))

        # Check if any alignments were found
        if blast_record.alignments:
            # The first alignment is the best hit
            top_hit = blast_record.alignments[0]
            print("\n--- Search Complete ---")
            print("The provided sequence most likely represents the following protein:")
            # The title of the hit contains the protein name, organism, and other info
            print(f"Description: {top_hit.title}")
            print(f"Length: {top_hit.length}")
            # The E-value (Expect value) describes the number of hits one can "expect"
            # to see by chance. The lower the E-value, the more significant the match.
            top_hsp = top_hit.hsps[0]
            print(f"E-value: {top_hsp.expect}")
            
            # Return the title for the final answer format
            return top_hit.title
        else:
            print("No significant matches found for the provided sequence.")
            return "No significant matches found."

    except Exception as e:
        print(f"\nAn error occurred during the search: {e}")
        print("This might be due to a network issue or a problem with the NCBI servers.")
        print("Please ensure you have an active internet connection and try again.")
        return "Error during search."

if __name__ == '__main__':
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
    
    protein_name = identify_protein_sequence(sequence)
    # The final answer format requires enclosing the result in <<<>>>
    # We will use the full description from the top BLAST hit.
    print(f"\n<<<{protein_name}>>>")
