# First, ensure you have biopython installed:
# pip install biopython

try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Biopython library not found.")
    print("Please install it using: pip install biopython")
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

# Remove whitespace and newlines to create a single sequence string
protein_sequence = "".join(sequence_data.split())

print(f"Identifying protein for sequence starting with: {protein_sequence[:30]}...")
print("Performing BLAST search against NCBI's 'nr' database. This may take a minute...")

try:
    # Use qblast to send the sequence to the NCBI BLAST server
    # program='blastp' for protein-protein comparison
    # database='nr' for the non-redundant protein database
    result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence)

    # Parse the XML results returned by the server
    # We use read() because we expect only one query result
    blast_record = NCBIXML.read(result_handle)

    # Check if any matching alignments were found
    if blast_record.alignments:
        # The first alignment is the best hit
        top_alignment = blast_record.alignments[0]
        
        # The title of the hit contains the protein's description
        protein_name = top_alignment.title
        
        print("\n--- Top Match Found ---")
        print(f"Description: {protein_name}")
        
    else:
        print("\nNo significant match found for the provided sequence.")

finally:
    # It's good practice to close the handle
    if 'result_handle' in locals() and result_handle:
        result_handle.close()

except Exception as e:
    print(f"\nAn error occurred during the BLAST search: {e}")
    print("This may be due to network issues or a problem with the NCBI server.")
