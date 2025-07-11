# The Biopython library is required. If not installed, run: pip install biopython
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Biopython library not found.")
    print("Please install it by running: pip install biopython")
    exit()

# The amino acid sequence provided by the user
sequence_data = """
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

# Clean up the sequence: remove spaces and newline characters
cleaned_sequence = sequence_data.replace(" ", "").replace("\n", "")

print("Identifying protein from sequence:")
print(cleaned_sequence[:60] + "...") # Print the start of the sequence
print("\nPerforming BLAST search against the NCBI 'nr' database. This may take a moment...")

try:
    # Use Biopython to send the request to the NCBI BLAST server
    # program='blastp' for protein vs protein database
    # database='nr' for non-redundant protein database
    # sequence=cleaned_sequence is our query
    result_handle = NCBIWWW.qblast("blastp", "nr", cleaned_sequence)

    # Parse the XML result returned by the server
    blast_record = NCBIXML.read(result_handle)

    # Check if any results were found
    if blast_record.alignments:
        print("\n--- Top Match Found ---")
        # The first alignment is the best hit
        top_alignment = blast_record.alignments[0]
        # The first hsp (high-scoring pair) contains the alignment details
        top_hsp = top_alignment.hsps[0]

        print(f"Description: {top_alignment.title}")
        print(f"Sequence Length: {top_alignment.length}")
        print(f"E-value (Expect): {top_hsp.expect}")
        print(f"Score: {top_hsp.score}")
        identity_percentage = 100 * top_hsp.identities / top_hsp.align_length
        print(f"Identity: {top_hsp.identities}/{top_hsp.align_length} ({identity_percentage:.2f}%)")
    else:
        print("\nNo significant match was found for the provided sequence.")

except Exception as e:
    print(f"\nAn error occurred: {e}")
    print("This could be due to a network issue or a problem with the NCBI services.")

finally:
    # Close the handle to release the connection
    if 'result_handle' in locals() and result_handle:
        result_handle.close()