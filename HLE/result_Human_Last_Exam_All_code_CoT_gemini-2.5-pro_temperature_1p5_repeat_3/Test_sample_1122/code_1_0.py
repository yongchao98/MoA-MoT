import sys
import os

# --- Step 1: Check for Biopython and install if not present ---
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("The 'biopython' library is required but not found.")
    # Attempt to install Biopython
    try:
        import pip
        print("Attempting to install biopython...")
        # Use a subprocess to ensure clean installation
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython"])
        print("\nInstallation successful! Please run the script again.")
    except Exception as e:
        print(f"\nFailed to auto-install biopython. Error: {e}")
        print("Please install it manually by running this command in your terminal:")
        print("pip install biopython")
    sys.exit()

# --- Step 2: Define the amino acid sequence ---
# The sequence is cleaned by removing spaces and newline characters.
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
protein_sequence = "".join(sequence_data.split())

# --- Step 3: Perform BLAST search ---
print("Identifying protein by searching the NCBI database using BLAST...")
print("This may take a moment, please be patient.")

try:
    # Use qblast to query the NCBI BLAST server.
    # program: 'blastp' (protein vs protein)
    # database: 'nr' (non-redundant protein database)
    # sequence: our protein sequence
    result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence)

    # --- Step 4: Parse the results ---
    # We sent one query, so we expect one BLAST record in the result.
    blast_record = NCBIXML.read(result_handle)

    # --- Step 5: Print the top hit information ---
    if blast_record.alignments:
        # The first alignment is the best match.
        top_hit = blast_record.alignments[0]
        print("\n--- Top Match Found ---")
        # The title contains the full description of the matched protein.
        print(f"Description: {top_hit.title}")
        print(f"Accession ID: {top_hit.accession}")
        # The E-value (Expect value) indicates the statistical significance of the match.
        # A very low E-value (close to 0) means the match is not by chance.
        print(f"E-value: {top_hit.hsps[0].expect}")
    else:
        print("\nNo significant match was found for the provided sequence.")

finally:
    # Always close the handle.
    if 'result_handle' in locals() and result_handle:
        result_handle.close()
