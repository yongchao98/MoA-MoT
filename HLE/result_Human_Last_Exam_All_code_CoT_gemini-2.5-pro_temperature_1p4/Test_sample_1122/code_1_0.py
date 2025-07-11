import sys

# Step 1: Check for Biopython and provide installation instructions if missing.
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Error: The 'biopython' library is required to run this script.")
    print("Please install it by running: pip install biopython")
    sys.exit(1)

# The amino acid sequence to identify
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
# Clean up the sequence by removing whitespace and newlines
sequence_string = "".join(sequence_data.strip().split())

print("Submitting BLASTP search to NCBI. This may take a minute...")

result_handle = None
try:
    # Step 2 & 3: Submit the blastp search to the NCBI server
    # - program='blastp' (protein-protein comparison)
    # - database='nr' (non-redundant protein database)
    # - sequence=the query sequence
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence_string)

    # Step 4: Parse the XML results returned by the server
    blast_record = NCBIXML.read(result_handle)

    # Step 5: Check if any matching alignments were found
    if blast_record.alignments:
        # The first alignment in the record is the best match
        top_alignment = blast_record.alignments[0]
        # The title of the alignment contains the description
        protein_description = top_alignment.title
        # The first HSP (High-scoring Segment Pair) contains the stats
        top_hsp = top_alignment.hsps[0]

        # Step 6: Print the results in a user-friendly format
        print("\n--- Search Complete: Top Match Found ---")
        
        # Extract the protein name and organism from the description string
        try:
            name_part = protein_description.split('|')[-1].split('[')[0].strip()
            organism_part = protein_description.split('[')[-1].split(']')[0].strip()
            print(f"Protein: {name_part}")
            print(f"Organism: {organism_part}")
        except IndexError:
            # Fallback for unexpected description formats
            print(f"Full Description: {protein_description}")
            
        print(f"\nSignificance (E-value): {top_hsp.expect}")
        print(f"Identity Match: {top_hsp.identities}/{top_hsp.align_length} ({100 * top_hsp.identities / top_hsp.align_length:.2f}%)")
    else:
        print("No significant matches were found for the provided sequence.")

except Exception as e:
    print(f"An error occurred during the BLAST search: {e}")
    print("Please check your internet connection and try again.")

finally:
    # Ensure the result handle is closed to release the connection
    if result_handle:
        result_handle.close()