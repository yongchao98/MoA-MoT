import io
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_protein_sequence():
    """
    Identifies a protein by performing a BLAST search with its amino acid sequence.
    """
    # The amino acid sequence provided by the user
    protein_sequence = (
        "MAEQVALSRT QVCGILREEL FQGDAFHQSD THIFIIMGAS GDLAKKKIYP "
        "TIWWLFRDGL LPENTFIVGY ARSRLTVADI RKQSEPFFKA TPEEKLKLED "
        "FFARNSYVAG QYDDAASYQR LNSHMNALHL GSQANRLFYL ALPPTVYEAV "
        "TKNIHESCMQ IGWNRIIVEK PFGRDLQSSD RLSNHISSLF REDQIYRIDH "
        "YLGKEMVQNL MVLRFANRIF GPIWNRDNIA CVILTFKEPF GTEGRGGYFD "
        "EFGIIRDVMQ NHLLQMLCLV AMEKPASTNS DDVRDEKVKV LKCISEVQAN "
        "NVVLGQYVGN PDGEGEATKG YLDDPTVPRG STTATFAAVV LYVENERWDG "
        "VPFILRCGKA LNERKAEVRL QFHDVAGDIF HQQCKRNELV IRVQPNEAVY "
        "TKMMTKKPGM FFNPEESELD LTYGNRYKNV KLPDAYERLI LDVFCGSQMH "
        "FVRSDELREA WRIFTPLLHQ IELEKPKPIP YIYGSRGPTE ADELMKRVGF "
        "QYEGTYKWVN PHKL"
    )

    # Clean up the sequence by removing spaces and newlines
    cleaned_sequence = "".join(protein_sequence.split())

    print("Performing BLAST search against the NCBI database. This may take a moment...")

    try:
        # Perform the BLAST search using the blastp program against the nr database
        result_handle = NCBIWWW.qblast("blastp", "nr", cleaned_sequence)

        # The result is in XML format. We read it and parse it.
        blast_xml_string = result_handle.read()
        result_handle.close()

        # Use io.StringIO to treat the XML string as a file for parsing
        blast_record = NCBIXML.read(io.StringIO(blast_xml_string))

        # Check if any alignments were found
        if blast_record.alignments:
            # Get the first (best) alignment
            top_alignment = blast_record.alignments[0]
            protein_title = top_alignment.title
            top_hsp = top_alignment.hsps[0]
            e_value = top_hsp.expect

            print("\n--- Top BLAST Hit ---")
            # The title contains the description of the matched protein
            print(f"Description: {protein_title}")
            print(f"E-value: {e_value}")
            print("\nThe E-value (Expect value) represents the number of hits one can expect to see by chance. A very low E-value (close to zero) indicates a highly significant match.")
        else:
            print("No significant protein matches were found for the given sequence.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an internet connection and the 'biopython' library is installed (`pip install biopython`).")

if __name__ == "__main__":
    identify_protein_sequence()