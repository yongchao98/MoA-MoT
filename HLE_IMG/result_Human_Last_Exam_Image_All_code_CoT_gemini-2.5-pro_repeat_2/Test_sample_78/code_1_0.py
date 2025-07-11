import re

def analyze_genotypes():
    """
    Analyzes the provided genetic data and gel image to determine the number of
    progeny for each genotype.
    """
    # Step 1: Define the initial biological data
    wt_orf = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    f_primer = "TTTTACGCGCTCTTCGTTTT"
    # Reverse complement of the reverse primer TTTTCCCTTGTCCACGAAAC
    r_primer_rc = "GTTTCGTGGACAAGGGAAAA"
    sfa_ni_site = "GCATC"

    # Step 2: In-silico PCR to find the amplicon
    print("--- Step 1: In-silico Analysis ---")
    f_primer_start = wt_orf.find(f_primer)
    r_primer_end = wt_orf.find(r_primer_rc) + len(r_primer_rc)
    wt_amplicon = wt_orf[f_primer_start:r_primer_end]
    amplicon_size = len(wt_amplicon)

    print(f"The PCR primers amplify a {amplicon_size} bp fragment from the wild-type gene.")

    # Step 3: Analyze restriction sites
    if sfa_ni_site in wt_amplicon or sfa_ni_site[::-1].translate(str.maketrans("ATCG", "TAGC")) in wt_amplicon:
        print("The wild-type amplicon contains an SfaNI site.")
    else:
        print("The wild-type amplicon does NOT contain an SfaNI site.")
        print("This means the C-to-A mutation must CREATE an SfaNI site for this assay to work.")

    # Step 4: Predict gel patterns based on the RFLP logic
    print("\n--- Step 2: Predicting Gel Patterns ---")
    print("Based on the principle of this RFLP (Restriction Fragment Length Polymorphism) assay:")
    print(f"  - Homozygous Wild-Type (+/+) will have one uncut band of {amplicon_size} bp.")
    print("  - Homozygous Mutant (-/-) will have its PCR product cut by SfaNI, resulting in smaller bands.")
    print("  - Heterozygote (+/-) will show both the uncut {amplicon_size} bp band and the smaller cut bands.")

    # Step 5: Interpret the provided gel image
    print("\n--- Step 3: Interpreting the Gel Image ---")
    print("The gel shows 16 sample lanes. The ladder on the right has its lowest band at 250 bp.")
    print(f"The uncut {amplicon_size} bp band corresponds to the high band seen just below the 250 bp marker.")
    print("The smaller, cut bands correspond to the low band(s) seen on the gel.")
    
    # Counts are determined by visual inspection of the gel image provided.
    # Pattern 1: Single high band (~238 bp) -> Homozygous Wild-Type. Lanes 4, 15.
    homozygous_wt_count = 2
    # Pattern 2: Both high and low bands -> Heterozygote. Lanes 1, 2, 3, 5, 8, 9, 10, 12, 13, 16.
    heterozygote_count = 10
    # Pattern 3: Only low band(s) -> Homozygous Mutant. Lanes 6, 7, 11, 14.
    homozygous_mutant_count = 4

    print("\nBy counting the lanes with each pattern, we find:")
    print(f"  - {homozygous_wt_count} lanes show the Homozygous Wild-Type pattern.")
    print(f"  - {heterozygote_count} lanes show the Heterozygous pattern.")
    print(f"  - {homozygous_mutant_count} lanes show the Homozygous Mutant pattern.")

    # Step 6: Final Answer
    print("\n--- Final Answer ---")
    print("The number of homozygous wild-type larvae (A) is: 2")
    print("The number of heterozygous larvae (B) is: 10")
    print("The number of homozygous mutant larvae (C) is: 4")
    print(f"\nThe final result in A/B/C format is: {homozygous_wt_count}/{heterozygote_count}/{homozygous_mutant_count}")

analyze_genotypes()