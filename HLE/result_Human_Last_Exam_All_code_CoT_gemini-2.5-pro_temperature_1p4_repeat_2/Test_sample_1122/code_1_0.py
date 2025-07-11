import sys

# Step 1: Check if BioPython is installed
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Error: The 'biopython' library is required to run this script.")
    print("Please install it by running: pip install biopython")
    sys.exit(1)

# Step 2: The amino acid sequence to identify
sequence_data = "MAEQVALSRTQVCGILREELFQGDAFHQSDTHIFIIMGASGDLAKKKIYPTIWWLFRDGLLPENTFIVGYARSRLTVADIRKQSEPFFKATPEEKLKLEDFFARNSYVAGQYDDAASYQRLNSHMNALHLGSQANRLFYLALPPTVYEAVTKNIHESCMQIGWNRIIVEKPFGRDLQSSDRLSNHISSLFREDQIYRIDHYLGKEMVQNLMVLRFANRIFGPIWNRDNIACVILTFKEPFGTEGRGGYFDEFGIIRDVMQNHLLQMLCLVAMEKPASTNSDDVRDEKVKVLKCISEVQANNVVLGQYVGNPDGEGEATKGYLDDPTVPRGSTATTFAAVVLYVENERWDGVPFILRCGKALNERKAEVRLQFHDVAGDIFHQQCKRNELVIRVQPNEAVYTKKMMTKKPGMFFNPEESELDLTYGNRYKNVKLPDAYERLILDVFCGSQMHFVRSDELREAWRIFTPLLHQIELEKPKPIPYIYGSRGPTEADELMKRVGFQYEGTYKWVNPHKL"

print("Identifying protein sequence...")
print("Running BLAST search against NCBI's database. This may take a moment.")

try:
    # Step 3: Perform the BLAST search using NCBI's online service
    # We use 'blastp' for protein-protein comparison and 'nr' for the non-redundant database
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence_data)

    # Step 4: Parse the results from the XML format returned by the server
    blast_record = NCBIXML.read(result_handle)

    # Step 5: Check if any matches were found
    if blast_record.alignments:
        # The first alignment is the best match
        top_alignment = blast_record.alignments[0]
        
        # The title of the alignment contains the protein's description
        protein_title = top_alignment.title
        
        print("\n--- Protein Identification Result ---")
        print("The most likely match for the provided sequence is:")
        print(protein_title)
        
        # Also show some alignment statistics for context
        top_hsp = top_alignment.hsps[0] # High-scoring Segment Pair
        print(f"\nMatch Statistics:")
        print(f"  Score: {top_hsp.score}")
        print(f"  E-value (Expect): {top_hsp.expect}") # Lower is better
        print(f"  Identities: {top_hsp.identities}/{top_hsp.align_length} ({(top_hsp.identities/top_hsp.align_length*100):.2f}%)")

    else:
        print("\nNo significant similarity was found for the provided sequence.")

except Exception as e:
    print(f"\nAn error occurred during the BLAST search: {e}")
    print("This could be a temporary network issue or a problem with the NCBI server.")

finally:
    if 'result_handle' in locals() and result_handle:
        result_handle.close()
