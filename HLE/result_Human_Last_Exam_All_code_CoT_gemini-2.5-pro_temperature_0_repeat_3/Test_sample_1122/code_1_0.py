import io
from Bio.Blast import NCBIWWW, NCBIXML

def identify_protein_sequence(sequence_data):
    """
    Identifies a protein sequence using NCBI's BLAST service.

    Args:
        sequence_data (str): A string containing the amino acid sequence.
    """
    # Clean the sequence by removing whitespace and newlines
    cleaned_sequence = "".join(sequence_data.split())

    print("Identifying protein using NCBI BLAST...")
    print("This may take a few moments as it queries a remote database.")

    try:
        # Perform the BLAST search against the non-redundant (nr) protein database
        result_handle = NCBIWWW.qblast("blastp", "nr", cleaned_sequence)

        # Read the XML result from the handle
        blast_xml_string = result_handle.read()
        result_handle.close()

        # To avoid re-running the query, we parse the saved XML string
        # Using io.StringIO to treat the string as a file
        blast_record = NCBIXML.read(io.StringIO(blast_xml_string))

        # Check if any alignments were found
        if blast_record.alignments:
            # The first alignment is the best match
            top_alignment = blast_record.alignments[0]
            top_hsp = top_alignment.hsps[0]

            print("\n--- Top Match Found ---")
            # The title contains the protein name, organism, and accession number
            print(f"Description: {top_alignment.title}")
            print(f"Sequence Length: {top_alignment.length}")
            # The E-value (Expect value) indicates the statistical significance of the match
            print(f"E-value: {top_hsp.expect}")
            print(f"Identity: {top_hsp.identities}/{top_hsp.align_length} ({100 * top_hsp.identities / top_hsp.align_length:.2f}%)")
            print("-----------------------")
            print("\nAn E-value close to 0.0 indicates a highly significant match, meaning it is not due to random chance.")

        else:
            print("No significant matches were found for the provided sequence.")

    except Exception as e:
        print(f"\nAn error occurred: {e}")
        print("This could be due to a network issue or a problem with the NCBI BLAST service.")
        print("Please ensure you have the 'biopython' library installed (`pip install biopython`) and an active internet connection.")

if __name__ == '__main__':
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