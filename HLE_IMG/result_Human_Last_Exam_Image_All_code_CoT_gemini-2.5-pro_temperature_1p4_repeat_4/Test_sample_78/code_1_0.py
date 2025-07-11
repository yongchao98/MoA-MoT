import re

def solve_genotyping_puzzle():
    """
    This script performs an in silico analysis to determine the expected RFLP patterns
    and uses this to interpret the provided gel image to genotype the zebrafish progeny.
    """

    # 1. Define Input Sequences and Enzyme
    wt_orf = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    rev_primer = "TTTTCCCTTGTCCACGAAAC"
    sfaNI_site = "GCATC"

    print("Step 1: In Silico PCR Analysis")

    # Find primer binding sites
    fwd_primer_start = wt_orf.find(fwd_primer)
    
    # Create reverse complement of the reverse primer
    rev_primer_rc = rev_primer.translate(str.maketrans("ATGC", "TACG"))[::-1]
    rev_primer_end = wt_orf.find(rev_primer_rc) + len(rev_primer_rc)

    # Extract PCR product
    pcr_product_wt = wt_orf[fwd_primer_start:rev_primer_end]
    pcr_product_size = len(pcr_product_wt)

    print(f"The forward primer binds at position {fwd_primer_start + 1}.")
    print(f"The reverse primer binds at position {rev_primer_end}.")
    print(f"The predicted full-length PCR product is {pcr_product_size} bp long.")
    print("-" * 30)

    print("Step 2 & 3: In Silico Restriction Digest Analysis (WT vs Mutant)")

    # Find SfaNI site in the WT PCR product
    match = re.search(sfaNI_site, pcr_product_wt)
    if match:
        cut_position = match.start() + len(sfaNI_site)
        fragment1_size = cut_position
        fragment2_size = pcr_product_size - cut_position
        print(f"The SfaNI recognition site '{sfaNI_site}' was found in the wild-type PCR product.")
        print("This implies the wild-type allele is digestible.")
        print(f"Digestion of the WT allele results in two fragments of sizes: {fragment1_size} bp and {fragment2_size} bp.")
    else:
        print("No SfaNI site found in the WT sequence. There is a contradiction in the problem statement.")
        return

    print("\nBased on the gel image showing three patterns, the C->A mutation must destroy this SfaNI site.")
    print(f"Therefore, the mutant allele is indigestible, resulting in a single uncut fragment of {pcr_product_size} bp.")
    print("-" * 30)
    
    print("Step 4: Predicted Genotype Patterns on the Gel")
    print(f"Homozygous Wild-Type (WT/WT) pattern: Two bands at {fragment1_size} bp and {fragment2_size} bp.")
    print(f"Homozygous Mutant (MUT/MUT) pattern: One band at {pcr_product_size} bp.")
    print(f"Heterozygous (WT/MUT) pattern: Three bands at {pcr_product_size} bp, {fragment1_size} bp, and {fragment2_size} bp.")
    print("-" * 30)
    
    print("Step 5: Analysis of the Gel Image")
    print("By comparing the predicted patterns with the bands on the gel (and using the ladder where the lowest band is 250 bp):")
    print(f"- The pattern with one large band (~{pcr_product_size} bp) corresponds to the homozygous mutant.")
    print(f"- The pattern with two small bands (~{fragment1_size} bp and ~{fragment2_size} bp) corresponds to the homozygous wild-type.")
    print(f"- The pattern with all three bands corresponds to the heterozygote.")
    
    # Counting lanes from the provided image
    wt_count = 5  # Lanes with 2 bands: 4, 5, 8, 14, 15
    het_count = 10 # Lanes with 3 bands: 1, 3, 6, 7, 9, 10, 11, 12, 13, 16
    mut_count = 1  # Lane with 1 band: 2
    
    print("\nCounting the lanes for each pattern on the gel gives:")
    print(f"Number of lanes with the WT/WT pattern: {wt_count}")
    print(f"Number of lanes with the WT/MUT (heterozygous) pattern: {het_count}")
    print(f"Number of lanes with the MUT/MUT pattern: {mut_count}")
    print("-" * 30)

    print("Step 6: Final Answer")
    print("The number of homozygous wild-type, heterozygous, and homozygous mutant larvae are:")
    print(f"Number of homozygous wild-type larvae (A): {wt_count}")
    print(f"Number of heterozygous larvae (B): {het_count}")
    print(f"Number of homozygous mutant larvae (C): {mut_count}")
    
    final_answer = f"{wt_count}/{het_count}/{mut_count}"
    print(f"\nFinal answer in A/B/C format: {final_answer}")


solve_genotyping_puzzle()