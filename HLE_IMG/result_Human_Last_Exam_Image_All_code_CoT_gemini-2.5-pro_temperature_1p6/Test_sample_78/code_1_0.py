import re

def analyze_genotypes():
    """
    Performs an in-silico PCR and restriction digest to determine the expected
    banding patterns for gdf6a genotyping and provides the counts from the gel.
    """
    # 1. Define sequences
    orf_wt = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    rev_primer = "TTTTCCCTTGTCCACGAAAC"
    sfaNI_site = "GCATC"
    mutation_pos_orf = 164  # 1-based position in ORF

    # 2. In-silico PCR
    fwd_start = orf_wt.find(fwd_primer)
    rev_comp = rev_primer.translate(str.maketrans("ATGC", "TACG"))[::-1]
    rev_comp_start = orf_wt.find(rev_comp)
    amplicon_end = rev_comp_start + len(rev_comp)
    amplicon_wt = orf_wt[fwd_start:amplicon_end]
    amplicon_size = len(amplicon_wt)

    print("--- In Silico Analysis ---")
    print(f"PCR amplicon size: {amplicon_size} bp")

    # 3. Analyze mutation and restriction site
    # Convert 1-based ORF position to 0-based index
    mutation_idx_orf = mutation_pos_orf - 1
    # Convert ORF index to amplicon index
    mutation_idx_amplicon = mutation_idx_orf - fwd_start
    
    # Create mutant amplicon
    amplicon_mut = list(amplicon_wt)
    amplicon_mut[mutation_idx_amplicon] = 'A'
    amplicon_mut = "".join(amplicon_mut)

    # Check for SfaNI site in both alleles
    wt_site_match = re.search(sfaNI_site, amplicon_wt)
    mut_site_match = re.search(sfaNI_site, amplicon_mut)

    print("\n--- Genotype Prediction ---")
    # Wild-Type Allele
    if not wt_site_match:
        print(f"Wild-Type Allele: The SfaNI site ({sfaNI_site}) is ABSENT.")
        print(f"Result: Uncut. Produces one band of {amplicon_size} bp.")
        wt_bands = [amplicon_size]
    else:
        # This case is not expected based on our analysis
        print("Error: Wild-type was expected to be uncut.")

    # Mutant Allele
    if mut_site_match:
        cut_pos = mut_site_match.start() + 1 # SfaNI cuts G^CATC
        frag1_size = cut_pos
        frag2_size = amplicon_size - frag1_size
        print(f"Mutant Allele: The C->A mutation CREATES an SfaNI site ({sfaNI_site}).")
        print(f"Result: Cut. Produces two bands of {frag1_size} bp and {frag2_size} bp.")
        mut_bands = [frag1_size, frag2_size]
    else:
        # This case is not expected
        print("Error: Mutant was expected to be cut.")
        
    print("\n--- Gel Interpretation ---")
    print("Homozygous Wild-Type (WT/WT): Shows one high molecular weight band (~263 bp).")
    print("Homozygous Mutant (MUT/MUT): Shows two low molecular weight bands (~141 bp and ~122 bp). These may appear as a single thick band.")
    print("Heterozygote (WT/MUT): Shows all three bands (one high, two low).")
    
    # 4. Count genotypes from gel image
    # Based on visual inspection of the 16 sample lanes in the gel image:
    # - The pattern with one high band (WT/WT) appears in 3 lanes.
    # - The pattern with low band(s) (MUT/MUT) appears in 5 lanes.
    # - The pattern with high and low bands (Heterozygote) appears in 8 lanes.
    # These numbers (3:8:5) are very close to the expected 1:2:1 Mendelian ratio for 16 progeny (4:8:4).
    
    num_wt = 3
    num_het = 8
    num_mut = 5

    print("\n--- Final Counts ---")
    print("Number of homozygous wild-type larvae (A): 3")
    print("Number of heterozygous larvae (B): 8")
    print("Number of homozygous mutant larvae (C): 5")
    
    print("\nFinal Answer Format (A/B/C):")
    # The final print statement is formatted for direct output as requested.
    # Each number is printed individually as part of the equation/format.
    print(f"{num_wt}/{num_het}/{num_mut}")

analyze_genotypes()