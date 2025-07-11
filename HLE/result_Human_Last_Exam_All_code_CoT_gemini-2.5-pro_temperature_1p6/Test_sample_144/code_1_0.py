import collections

def solve_genetics_sequences():
    """
    Calculates the number of unique chromosome sequences possible in the F3 generation.
    """
    # Number of SNPs on the autosome
    n_snps = 5

    # Step 1: Define the parental (P0) strains and their gametes.
    # Using 'A' and 'B' to represent the alleles from the two inbred strains.
    p0_strain_a = "A" * n_snps
    p0_strain_b = "B" * n_snps
    s_p0_gametes = {p0_strain_a, p0_strain_b}
    num_p0_gametes = len(s_p0_gametes)
    print(f"The two parental (P0) sequences are: {p0_strain_a} and {p0_strain_b}")
    print(f"Number of unique P0 gametes: {num_p0_gametes}")
    print("-" * 20)

    # Helper function to generate recombinant sequences from a single crossover.
    def get_recombinants(seq1, seq2):
        recombinant_seqs = set()
        n = len(seq1)
        # A single crossover can occur at n-1 locations between the SNPs.
        for i in range(1, n):
            prefix1, suffix1 = seq1[:i], seq1[i:]
            prefix2, suffix2 = seq2[:i], seq2[i:]
            recombinant_seqs.add(prefix1 + suffix2)
            recombinant_seqs.add(prefix2 + suffix1)
        return recombinant_seqs

    # Step 2: Calculate the unique gametes produced by the F1 generation.
    # These gametes form the F2 generation.
    s_f1_gametes = s_p0_gametes.copy()
    p0_list = list(s_p0_gametes)
    f1_recombinants = get_recombinants(p0_list[0], p0_list[1])
    s_f1_gametes.update(f1_recombinants)
    
    num_crossover_points = n_snps - 1
    num_recombinant_F1 = len(f1_recombinants)
    num_f1_gametes = len(s_f1_gametes)
    
    print(f"F1 individuals are heterozygous ({p0_strain_a}/{p0_strain_b}).")
    print("They produce gametes with one crossover.")
    print(f"Number of crossover points between {n_snps} SNPs = {n_snps} - 1 = {num_crossover_points}")
    print(f"Number of unique F1 gametes = (Parental) + (Recombinant)")
    print(f"Number of unique F1 gametes = {num_p0_gametes} + 2 * {num_crossover_points} = {num_f1_gametes}")
    print(f"Set of F1 gametes: {sorted(list(s_f1_gametes))}")
    print("-" * 20)

    # Step 3: Calculate the unique gametes produced by the F2 generation.
    # These sequences will be found in the F3 individuals.
    s_f2_gametes = s_f1_gametes.copy()
    f1_gametes_list = list(s_f1_gametes)

    # Iterate through all pairs of F1 gametes to find all possible recombinants.
    for i in range(len(f1_gametes_list)):
        for j in range(i, len(f1_gametes_list)):
            f2_recombinants = get_recombinants(f1_gametes_list[i], f1_gametes_list[j])
            s_f2_gametes.update(f2_recombinants)

    final_count = len(s_f2_gametes)

    # The final number should be 2^n_snps, as any combination becomes possible.
    print("F2 individuals produce gametes by recombining any pair of the F1 gamete types.")
    print("This process generates all possible combinations of the original SNPs.")
    final_equation_part_1 = "2"
    final_equation_part_2 = str(n_snps)
    print(f"Final number of unique sequences = {final_equation_part_1}^{final_equation_part_2} = {final_count}")
    
solve_genetics_sequences()
