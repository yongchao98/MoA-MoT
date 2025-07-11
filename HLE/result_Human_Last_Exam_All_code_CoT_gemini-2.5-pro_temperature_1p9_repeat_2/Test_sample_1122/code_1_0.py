import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from urllib.error import URLError

def identify_protein_sequence(sequence):
    """
    Identifies a protein sequence using NCBI BLAST.
    """
    print("Performing BLAST search. This may take a moment...")
    try:
        # Use qblast to submit the search to NCBI
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence, hitlist_size=1)
        
        # Read the single BLAST record from the result handle
        # NCBIXML.read() is used when we expect one and only one query result
        blast_record = NCBIXML.read(result_handle)
        
        # Check if any alignments were found
        if blast_record.alignments:
            # The top hit is the first alignment in the list
            top_hit = blast_record.alignments[0]
            print("\n--- Top BLAST Hit ---")
            print("Title:", top_hit.title)
            print("Length:", top_hit.length)
            # The e-value (expect value) describes the number of hits one can "expect"
            # to see by chance. The smaller the e-value, the more significant the match.
            top_hsp = top_hit.hsps[0] # High-scoring Segment Pair
            print("E-value:", top_hsp.expect)
            
            # Extract the common protein name for the final answer
            # The title often looks like 'gi|...|ref|NP_...| phosphoglycerate mutase 2 [Homo sapiens]'
            # We will try to extract the protein name and species.
            protein_name = top_hit.title.split('|')[-1].strip()
            # The answer format requires just the name, let's extract that.
            final_answer_name = protein_name.split('[')[0].strip()
            return final_answer_name

        else:
            print("No significant matches found for the sequence.")
            return "Unknown Protein"
            
    except URLError as e:
        print(f"Network error: Could not connect to NCBI. Please check your internet connection.", file=sys.stderr)
        print(f"Error details: {e}", file=sys.stderr)
        return "Error during search"
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        return "Error during search"

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
).replace(" ", "").replace("\n", "")

# Run the identification function
final_protein_name = identify_protein_sequence(protein_sequence)
# The final answer format is specified to be separate.
# So we'll print a concluding remark from the code itself.
print(f"\n<<<The provided sequence most likely represents: {final_protein_name}>>>")
