# Note: This script requires the Biopython library.
# You can install it using pip: pip install biopython

import sys

try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Error: Biopython library not found.")
    print("Please install it by running: pip install biopython")
    sys.exit(1)

# The amino acid sequence to be identified
sequence_data = (
    "MAEQVALSRTQVCGILREELFQGDAFHQSDTHIFIIMGASGDLAKKKIYP"
    "TIWWLFRDGLLPENTFIVGYARSRLTVADIRKQSEPFFKATPEEKLKLED"
    "FFARNSYVAGQYDDAASYQRLNSHMNALHLGSQANRLFYLALPPTVYEAV"
    "TKNIHESCMQIGWNRIIVEKPFGRDLQSSDRLSNHISSLFREDQIYRIDH"
    "YLGKEMVQNLLVLRFANRIFGPIWNRDNIACVILTFKEPFGTEGRGGYFD"
    "EFGIIRDVMQNHLLQMLCLVAMEKPASTNSDDVRDEKVKVLKCISEVQAN"
    "NVVLGQYVGNPDGEGEATKGYLDDPTVPRGSTTATFAAVVLYVENERWDG"
    "VPFILRCGKALNERKAEVRLQFHDVAGDIFHQQCKRNELVIRVQPNEAVY"
    "TKMMTKKPGMFFNPEESELDLTYGNRYKNVKLPDAYERLIILDVFCGSQMH"
    "FVRSDELREAWRIFTPLLHQIELEKPKPIPYIYGSRGPTEADELMKRVGF"
    "QYEGTYKWVNPHKL"
)

print("Performing protein identification using NCBI BLAST.")
print("This may take a minute or two...")

try:
    # Use the qblast function to submit the search to NCBI
    # 'blastp' is for protein sequence comparison
    # 'nr' is the non-redundant protein sequence database
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence_data)

    # The result is in XML format, which we need to parse
    blast_record = NCBIXML.read(result_handle)

    # Check if the search returned any results (alignments)
    if blast_record.alignments:
        # The results are pre-sorted, so the first alignment is the best one
        top_hit = blast_record.alignments[0]
        # The High-scoring Segment Pair (HSP) contains the match statistics
        top_hsp = top_hit.hsps[0]

        print("\n--- Top Match Found ---")
        # The title of the alignment contains the protein name and organism
        print(f"Description: {top_hit.title}")
        print(f"Sequence Length: {top_hit.length}")
        # The E-value (Expect value) tells us the significance of the match.
        # A value close to zero means the match is highly significant.
        print(f"E-value: {top_hsp.expect}")

        # Extract a more concise name for the final answer
        protein_info = "Unknown Protein"
        # The title often follows the format "description [Organism]"
        if '[' in top_hit.title and ']' in top_hit.title:
            protein_info = top_hit.title.split(' ', 1)[1] # Remove accession codes

        print(f"\nConclusion: The sequence represents '{protein_info}'.")

    else:
        print("\nNo significant match was found in the database for the given sequence.")

except Exception as e:
    print(f"\nAn error occurred: {e}")
    print("Please check your internet connection or try again later.")
finally:
    # It's good practice to close the handle
    if 'result_handle' in locals() and result_handle:
        result_handle.close()