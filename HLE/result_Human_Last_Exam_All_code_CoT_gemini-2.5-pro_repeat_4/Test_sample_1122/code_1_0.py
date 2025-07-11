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

# Clean up the sequence string (remove spaces and newlines)
cleaned_sequence = "".join(sequence_data.split())

print("Performing BLAST search... This may take a few moments.")

try:
    # Use Biopython to send the sequence to the NCBI BLAST server
    # 'blastp' is for protein-protein search
    # 'nr' is the non-redundant protein database
    result_handle = NCBIWWW.qblast("blastp", "nr", cleaned_sequence)

    # Parse the XML result returned by the server
    blast_record = NCBIXML.read(result_handle)

    # Check if any results were found
    if blast_record.alignments:
        # Get the first (and best) alignment
        top_alignment = blast_record.alignments[0]
        # Get the first High-scoring Segment Pair (HSP) for that alignment
        top_hsp = top_alignment.hsps[0]

        print("\n--- Top Match Found ---")
        # The title contains the description of the matched protein
        print(f"Protein Name: {top_alignment.title}")
        print(f"Sequence Length: {top_alignment.length}")
        # The E-value (Expect value) represents the probability of a match occurring by chance.
        # A very low E-value (close to 0) indicates a highly significant match.
        print(f"E-value: {top_hsp.expect}")
        print(f"Score: {top_hsp.score}")
    else:
        print("No matching protein was found for the given sequence.")

except Exception as e:
    print(f"An error occurred: {e}", file=sys.stderr)
    print("This may be due to a network issue or a problem with the NCBI server.", file=sys.stderr)

finally:
    if 'result_handle' in locals() and result_handle:
        result_handle.close()
