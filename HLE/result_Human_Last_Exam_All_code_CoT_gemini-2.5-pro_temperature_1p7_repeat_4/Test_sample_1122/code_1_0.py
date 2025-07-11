import sys

# It is a good practice to check for necessary libraries.
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Error: The 'BioPython' library is required to run this script.", file=sys.stderr)
    print("Please install it using: pip install biopython", file=sys.stderr)
    sys.exit(1)

def identify_protein_from_sequence(sequence_str):
    """
    Takes an amino acid sequence, performs a BLAST search, and prints the top hit.
    """
    # First, clean the sequence by removing any spaces or newline characters.
    sequence = "".join(sequence_str.strip().split())

    print("Submitting sequence to NCBI BLAST...")
    print("This search may take a moment to complete.")

    result_handle = None
    try:
        # Perform the BLASTp search against the 'nr' (non-redundant) protein database.
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence)

        # The result from BLAST is in XML format. We need to parse it.
        blast_record = NCBIXML.read(result_handle)

        # Check if the search returned any matching alignments.
        if blast_record.alignments:
            # The first alignment is the best match with the highest score.
            top_alignment = blast_record.alignments[0]
            # The 'hsps' (High-scoring Segment Pairs) contain scoring details.
            top_hsp = top_alignment.hsps[0]

            print("\n--- BLAST Search Complete ---")
            print(f"Top match found for the sequence.")
            # The title of the alignment contains the protein's name, source organism, and accession codes.
            print(f"Description: {top_alignment.title}")
            # The E-value (Expect value) represents the statistical significance of the match.
            # A very low E-value (close to zero) means the match is highly significant.
            print(f"E-value: {top_hsp.expect}")
        else:
            print("No significant matches were found for the provided sequence.")

    except Exception as e:
        print(f"\nAn error occurred: {e}", file=sys.stderr)
        print("Please check your internet connection and ensure NCBI services are available.", file=sys.stderr)
    finally:
        # It's important to close the result handle after we're done.
        if result_handle is not None:
            result_handle.close()

# The amino acid sequence provided in the prompt.
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

# Execute the function with the sequence.
if __name__ == '__main__':
    identify_protein_from_sequence(amino_acid_sequence)