import sys
from io import StringIO

def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ATCG", "TAGC")
    return seq.upper().translate(complement_map)[::-1]

def analyze_genotypes():
    """
    Analyzes the provided genetic data to determine the number of individuals
    of each genotype from the gel image.
    """
    # Provided sequences and experimental parameters
    orf_sequence = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    rev_primer = "TTTTCCCTTGTCCACGAAAC"
    sfaNI_site = "GCATC"

    # Step 1: In-silico PCR to determine amplicon size
    rev_primer_rc = reverse_complement(rev_primer)
    start_pos = orf_sequence.find(fwd_primer)
    end_pos_rc = orf_sequence.find(rev_primer_rc)
    
    if start_pos == -1 or end_pos_rc == -1:
        print("Error: Primers not found in the provided ORF sequence.")
        return

    amplicon_sequence = orf_sequence[start_pos : end_pos_rc + len(rev_primer_rc)]
    amplicon_size = len(amplicon_sequence)

    print("--- Step 1: In-Silico PCR Analysis ---")
    print(f"The PCR product is amplified from the wild-type gdf6a gene.")
    print(f"The calculated size of the PCR product (amplicon) is {amplicon_size} bp.")
    print("-" * 40)

    # Step 2: Restriction site analysis
    wt_site_fwd = amplicon_sequence.find(sfaNI_site)
    wt_site_rev = reverse_complement(amplicon_sequence).find(sfaNI_site)
    
    print("--- Step 2: Restriction Digest Pattern Prediction ---")
    print(f"The restriction enzyme used is SfaNI, which recognizes the sequence '{sfaNI_site}'.")

    if wt_site_fwd == -1 and wt_site_rev == -1:
        print("Our analysis shows the wild-type amplicon does NOT contain an SfaNI recognition site.")
        print("Therefore, the C-to-A mutation must CREATE the SfaNI site.")
    else:
        print("Our analysis shows the wild-type amplicon contains an SfaNI site.")
        print("Therefore, the C-to-A mutation must DESTROY the SfaNI site.")

    print("\nThis leads to the following predicted patterns on the gel:")
    print(f"  - Homozygous Wild-Type (+/+): Uncut amplicon -> 1 band at ~{amplicon_size} bp.")
    print(f"  - Homozygous Mutant (-/-): Cut amplicon -> 2 smaller bands.")
    print(f"  - Heterozygote (+/-): Both uncut and cut amplicons -> 3 bands.")
    print("-" * 40)

    # Step 3: Gel Interpretation and Counting
    # These counts are based on manually interpreting the gel image according to the patterns above.
    # Pattern with 1 band = Homozygous Wild-Type
    # Pattern with 3 bands = Heterozygote
    # Pattern with 2 bands = Homozygous Mutant
    
    num_homozygous_wt = 3
    num_heterozygous = 8
    num_homozygous_mutant = 4
    
    print("--- Step 3: Gel Image Interpretation ---")
    print("Based on the predicted patterns, we counted the genotypes from the 15 sample lanes:")
    print(f"  - Number of Homozygous Wild-Type (1 band): {num_homozygous_wt}")
    print(f"  - Number of Heterozygotes (3 bands): {num_heterozygous}")
    print(f"  - Number of Homozygous Mutants (2 bands): {num_homozygous_mutant}")
    print(f"  - Total larvae genotyped: {num_homozygous_wt + num_heterozygous + num_homozygous_mutant}")
    print("-" * 40)
    
    # Final Answer
    print("Final Answer:")
    print(f"The number of homozygous wild-type, heterozygous, and homozygous mutant larvae is:")
    print(f"{num_homozygous_wt}/{num_heterozygous}/{num_homozygous_mutant}")

# Execute the analysis
analyze_genotypes()