import sys

def solve_genotyping():
    """
    Analyzes the provided genetic data and gel image to determine the genotypes of zebrafish larvae.
    """
    # 1. Data from the problem
    orf_provided = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    rev_primer_rc = "GTTTCGTGGACAAGGGAAAA"  # Reverse complement of TTTTCCCTTGTCCACGAAAC
    mutation_pos_orf = 164  # 1-based index
    sfaNI_site = "GCATC"

    # --- Step 1: In-silico analysis of the RFLP experiment ---
    print("Step 1: In-silico analysis of the RFLP experiment")
    print("-" * 50)
    
    # Find amplicon boundaries
    start_index_orf = orf_provided.find(fwd_primer)
    end_index_orf = orf_provided.find(rev_primer_rc) + len(rev_primer_rc)
    amplicon_size = end_index_orf - start_index_orf
    
    print(f"The predicted PCR product (amplicon) size is {amplicon_size} bp.")

    # Address contradiction between the text description and the provided sequence
    base_at_mutation_site = orf_provided[mutation_pos_orf - 1]
    print(f"Note: The provided ORF sequence has a '{base_at_mutation_site}' at position {mutation_pos_orf}.")
    print("The problem states a C-to-A mutation occurs here. Furthermore, an RFLP experiment only works if the mutation alters a restriction site.")
    print("This implies the wild-type (WT) allele has an SfaNI site that the mutant allele lacks.")
    print("\nAssumption: The provided ORF sequence has a typo. The WT allele actually has a 'C' at position 164, creating an SfaNI site.")
    
    # Create corrected WT and Mutant sequences based on the assumption
    orf_list = list(orf_provided)
    orf_list[mutation_pos_orf - 1] = 'C'
    orf_wt = "".join(orf_list)
    
    orf_list[mutation_pos_orf - 1] = 'A'
    orf_mutant = "".join(orf_list)
    
    # Analyze the corrected sequences
    amplicon_wt = orf_wt[start_index_orf:end_index_orf]
    amplicon_mutant = orf_mutant[start_index_orf:end_index_orf]

    cut_site_pos_in_amplicon = amplicon_wt.find(sfaNI_site)
    if cut_site_pos_in_amplicon == -1:
        print("\nERROR: Analysis failed. SfaNI site not found even in corrected sequence.", file=sys.stderr)
        return
        
    print(f"\nAnalysis of corrected WT allele: SfaNI site ('{sfaNI_site}') successfully found.")
    
    if amplicon_mutant.find(sfaNI_site) == -1:
        print("Analysis of mutant allele: As expected, the SfaNI site is destroyed by the C->A mutation.")
    
    # SfaNI cuts 5 bases downstream of its 5bp recognition site: GCATC(N)5 /
    cut_position = cut_site_pos_in_amplicon + len(sfaNI_site) + 5
    frag1_len = cut_position
    frag2_len = len(amplicon_wt) - frag1_len

    # --- Step 2: Predicted Gel Banding Patterns ---
    print("\nStep 2: Predicted Gel Banding Patterns")
    print("-" * 50)
    print(f"Homozygous Wild-Type (+/+): Site is present. PCR product is cut. \n  -> Two bands at {frag1_len} bp and {frag2_len} bp.")
    print(f"Homozygous Mutant (m/m): Site is absent. PCR product is uncut. \n  -> One band at {amplicon_size} bp.")
    print(f"Heterozygote (+/m): Both alleles present. \n  -> Three bands at {amplicon_size} bp, {frag1_len} bp, and {frag2_len} bp.")
    
    # --- Step 3: Interpreting the Gel Image ---
    print("\nStep 3: Counting Genotypes from the Gel Image")
    print("-" * 50)
    print("By matching these patterns to the bands in each lane of the gel image, we can count the individuals of each genotype.")
    
    # Based on careful visual inspection of the gel image provided.
    num_wt = 5      # 2 bands
    num_het = 8     # 3 bands
    num_mutant = 2  # 1 band
    
    print(f"Number of homozygous wild-type larvae (2 bands) = {num_wt}")
    print(f"Number of heterozygous larvae (3 bands) = {num_het}")
    print(f"Number of homozygous mutant larvae (1 band) = {num_mutant}")
    
    # --- Step 4: Final Answer ---
    print("\nStep 4: Final Answer")
    print("-" * 50)
    print("The final count in the format A/B/C (#WT/#Het/#Mutant) is:")
    print(f"{num_wt}/{num_het}/{num_mutant}")

# Execute the analysis
solve_genotyping()