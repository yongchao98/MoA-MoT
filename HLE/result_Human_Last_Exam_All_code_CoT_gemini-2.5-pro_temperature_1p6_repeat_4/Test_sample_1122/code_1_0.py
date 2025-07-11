import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

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

# Clean up the sequence data by removing spaces and newlines
sequence_clean = "".join(sequence_data.split())

print("Identifying protein sequence...")
print("Please wait, this may take a minute as it contacts the NCBI BLAST server...")

try:
    # Run the BLAST search online
    # 'blastp' is for protein-protein comparison
    # 'nr' is the non-redundant protein database
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence_clean)
    
    # Parse the XML result
    blast_records = NCBIXML.parse(result_handle)
    
    # Get the first record (since we only submitted one query)
    blast_record = next(blast_records)
    
    # Check if any significant alignments were found
    if not blast_record.alignments:
        print("\nNo significant matches found for the provided sequence.")
    else:
        # The first alignment is the best hit
        top_alignment = blast_record.alignments[0]
        
        # Get the descriptive title of the top hit
        protein_title = top_alignment.title
        
        print("\n--- Top Match Found ---")
        print(f"The most likely protein is: {protein_title}")

except Exception as e:
    print(f"\nAn error occurred: {e}", file=sys.stderr)
    print("Please check your internet connection and ensure 'biopython' is installed (`pip install biopython`).", file=sys.stderr)

finally:
    # Clean up the handle if it was created
    if 'result_handle' in locals() and result_handle:
        result_handle.close()