import sys
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Biopython library not found.")
    print("Please install it using: pip install biopython")
    sys.exit(1)

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

# Clean up the sequence by removing spaces and newlines
cleaned_sequence = "".join(sequence_data.strip().split())

print("Performing protein BLAST search against the NCBI database...")
print("This may take a minute or two...")

try:
    # Perform the BLAST search using NCBI's web service
    # The 'nr' database is the non-redundant protein database
    result_handle = NCBIWWW.qblast(program="blastp", database="nr", sequence=cleaned_sequence)

    # Parse the XML result from BLAST
    blast_record = NCBIXML.read(result_handle)

    # Check if any significant alignments were found
    if not blast_record.alignments:
        print("No significant matches found for the provided sequence.")
    else:
        # Get the first alignment, which is the top hit
        top_alignment = blast_record.alignments[0]
        
        # The 'hit_def' attribute contains the description of the matched protein
        protein_description = top_alignment.hit_def
        e_value = top_alignment.hsps[0].expect

        print("\n--- Top Match Found ---")
        print(f"Protein: {protein_description}")
        print(f"E-value: {e_value:.2e} (a very low E-value indicates a highly significant match)")

except Exception as e:
    print(f"An error occurred during the BLAST search: {e}")
    print("This could be due to a network issue or a problem with the NCBI server.")

finally:
    # It's good practice to close the handle
    if 'result_handle' in locals() and result_handle:
        result_handle.close()
        