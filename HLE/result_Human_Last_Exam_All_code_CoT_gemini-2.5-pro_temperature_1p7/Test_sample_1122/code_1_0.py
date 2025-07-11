import sys
import io

# This script requires the biopython library.
# If you don't have it, you can install it by running: pip install biopython
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Error: The 'biopython' library is required to run this script.", file=sys.stderr)
    print("Please install it using: pip install biopython", file=sys.stderr)
    sys.exit(1)

# The amino acid sequence to be identified
sequence = "MAEQVALSRTQVCGILREELFQGDAFHQSDTHIFIIMGASGDLAKKKIYPTIWWLFRDGLLPENTFIVGYARSRLTVADIRKQSEPFFKATPEEKLKLEDFFARNSYVAGQYDDAASYQRLNSHMNALHLGSQANRLFYLALPPTVYEAVTKNIHESCMQIGWNRIIVEKPFGRDLQSSDRLSNHISSLFREDQIYRIDHYLGKEMVQNLMVLRFANRIFGPIWNRDNIACVILTFKEPFGTEGRGGYFDEFGIIRDVMQNHLLQMLCLVAMEKPASTNSDDVRDEKVKVLKCISEVQANNVVLGQYVGNPDGEGEATKGYLDDPTVPRGSTATFAAVVLYVENERWDGVPFILRCGKALNERKAEVRLQFHDVAGDIFHQQCKRNELVIRVQPNEAVYTKMMTKKPGMFFNPEESELDLTYGNRYKNVKLPDAYERLILDVFCGSQMHFVRSDELREAWRIFTPLLHQIELEKPKPIPIYGSRGPTEADELMKRVGFQYEGTYKWVNPHKL"

print("Querying the NCBI BLAST database to identify the protein...")
print("This process may take a minute or two, please wait.")

try:
    # Use NCBIWWW.qblast to perform the BLAST search
    # 'blastp' is for protein-protein BLAST
    # 'nr' is the non-redundant protein database
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence)

    # The result from NCBI is in XML format. We parse it directly.
    # Using io.StringIO allows parsing the results from a string in memory
    blast_records = NCBIXML.parse(result_handle)
    
    # We expect a single result record since we submitted one sequence
    blast_record = next(blast_records)

    # Check if any alignments were found
    if blast_record.alignments:
        # The first alignment object corresponds to the best (lowest E-value) hit
        top_alignment = blast_record.alignments[0]
        
        # The 'hit_def' attribute contains a clean definition of the protein and organism
        protein_definition = top_alignment.hit_def
        
        # The top High-scoring Segment Pair (HSP) contains statistics for the match
        top_hsp = top_alignment.hsps[0]
        
        print("\n--- Search Complete ---")
        print("Protein identified based on the highest-scoring match:")
        print(f"\nDescription: {protein_definition}")
        print(f"Match Score: {top_hsp.score}")
        print(f"E-value: {top_hsp.expect} (lower is better)")
        print(f"Identity: {top_hsp.identities}/{top_hsp.align_length} ({(100 * top_hsp.identities / top_hsp.align_length):.2f}%)")
        
        # Extract just the protein name, which is typically before the organism in brackets
        protein_name = protein_definition.split(' [')[0]
        
        print(f"\nBased on this, the provided sequence most likely represents the protein: {protein_name}")
        
    else:
        print("\nNo significant matches found for the provided sequence in the database.")

except Exception as e:
    print(f"\nAn error occurred during the BLAST search: {e}", file=sys.stderr)
    print("This may be due to a network issue or a temporary problem with the NCBI server.", file=sys.stderr)
