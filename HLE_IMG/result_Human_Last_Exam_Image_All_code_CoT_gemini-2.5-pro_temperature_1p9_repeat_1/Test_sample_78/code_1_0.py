import re

def solve_genotyping():
    """
    Analyzes zebrafish genotyping data to determine the number of homozygous wild-type,
    heterozygous, and homozygous mutant progeny.
    """
    # 1. Define sequences and primers
    wt_orf = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    rev_primer_rc = "GTTTCGTGGACAAGGGAAAA" # Reverse complement of TTTTCCCTTGTCCACGAAAC

    # 2. Locate amplicon in the ORF
    fwd_start = wt_orf.find(fwd_primer)
    rev_start = wt_orf.find(rev_primer_rc)
    pcr_end = rev_start + len(rev_primer_rc)
    
    wt_amplicon = wt_orf[fwd_start:pcr_end]
    pcr_product_size = len(wt_amplicon)

    # 3. Model the mutation
    # The mutation is a C->A change at ORF position 164 (string index 163).
    # First, we must find this position relative to the start of our amplicon.
    mutation_pos_in_orf = 163
    mutation_pos_in_amplicon = mutation_pos_in_orf - fwd_start
    
    # Create the mutant amplicon
    mutant_amplicon = list(wt_amplicon)
    mutant_amplicon[mutation_pos_in_amplicon] = 'A'
    mutant_amplicon = "".join(mutant_amplicon)

    # 4. Perform in silico restriction digest with SfaNI (GCATC)
    enzyme_site = "GCATC"
    wt_site = wt_amplicon.find(enzyme_site)
    mut_site = mutant_amplicon.find(enzyme_site)

    print("--- Step-by-Step Analysis ---")
    print(f"1. The PCR product size is {pcr_product_size} bp.")
    print(f"2. The mutation is a C->A change at position 164 of the ORF.")
    print(f"   Wild-type sequence at site: ...{wt_amplicon[mutation_pos_in_amplicon-5:mutation_pos_in_amplicon+5]}...")
    print(f"   Mutant sequence at site:  ...{mutant_amplicon[mutation_pos_in_amplicon-5:mutation_pos_in_amplicon+5]}...")

    if wt_site == -1 and mut_site != -1:
        print(f"3. The SfaNI recognition site '{enzyme_site}' is ABSENT in the wild-type allele.")
        print(f"   The mutation CREATES the SfaNI site in the mutant allele.")
        
        # Calculate fragment sizes. SfaNI cuts at GCA TC(N)5 /^...
        # The cut is 5 bases downstream from the end of the site.
        # This is a simplification; the exact cut pattern is GCATC(5/9),
        # but for fragment size calculation this is sufficient.
        cut_pos = mut_site + len(enzyme_site) + 4 # approximate cut position
        fragment1_size = cut_pos
        fragment2_size = pcr_product_size - fragment1_size
        
        print("\n--- Expected Gel Patterns ---")
        print(f"A) Homozygous Wild-Type (WT/WT): Uncut. One band at ~{pcr_product_size} bp.")
        print(f"B) Homozygous Mutant (Mut/Mut): Cut. Two bands at ~{fragment1_size} bp and ~{fragment2_size} bp.")
        print(f"C) Heterozygous (WT/Mut): Both uncut and cut DNA. Three bands at ~{pcr_product_size} bp, ~{fragment1_size} bp, and ~{fragment2_size} bp.")

    else:
        # This case handles other possibilities, though not expected here.
        print("Analysis Error: The mutation does not create or destroy the site as expected.")
        return

    # 5. Interpret the gel based on the derived patterns
    # These numbers come from visually inspecting the gel lanes.
    num_hom_wt = 3  # Lanes with one high band (e.g., lane 1, 14, 16)
    num_het = 8     # Lanes with three bands (e.g., lane 2, 5, 6, 7, 10-13)
    num_hom_mut = 4 # Lanes with two low bands (e.g., lane 3, 4, 8, 15)

    print("\n--- Final Counts from Gel Image ---")
    print(f"By observing the gel image and counting the lanes corresponding to each pattern:")
    print(f"Number of homozygous wild-type larvae (A) = {num_hom_wt}")
    print(f"Number of heterozygous larvae (B) = {num_het}")
    print(f"Number of homozygous mutant larvae (C) = {num_hom_mut}")
    
    print("\nFinal Answer Format:")
    print(f"{num_hom_wt}/{num_het}/{num_hom_mut}")


# Run the analysis
solve_genotyping()