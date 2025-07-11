import collections

def recombine(h1, h2):
    """
    Generates all possible single-crossover products between two haplotypes.
    A crossover can happen in the 4 intervals between the 5 SNPs.
    """
    products = set()
    # There are 5 SNPs, so 4 possible crossover locations (1, 2, 3, 4).
    # A crossover at location 'i' occurs after the i-th SNP.
    for i in range(1, 5):
        # Product 1: prefix from h1, suffix from h2
        prod1 = h1[:i] + h2[i:]
        products.add(prod1)
        # Product 2: prefix from h2, suffix from h1
        prod2 = h2[:i] + h1[i:]
        products.add(prod2)
    return products

def count_transitions(haplotype):
    """Counts the number of allele transitions (e.g., 0 to 1) in a sequence."""
    transitions = 0
    for i in range(len(haplotype) - 1):
        if haplotype[i] != haplotype[i+1]:
            transitions += 1
    return transitions

def solve_genetics_problem():
    """
    Calculates the number of unique autosome sequences in the F3 generation.
    """
    # Step 1: Define parental haplotypes
    parent_A = '00000'
    parent_B = '11111'

    # Step 2: Generate F1 gametes. According to the problem, every gamete
    # is a product of exactly one recombination.
    f1_gametes = recombine(parent_A, parent_B)

    # Step 3: Generate F2 gametes (these are the haplotypes found in the F3 generation).
    # F2 individuals are formed from pairs of F1 gametes. We simulate recombination
    # for every possible F2 genotype.
    f3_haplotypes = set()
    f1_gamete_list = list(f1_gametes)

    for i in range(len(f1_gamete_list)):
        for j in range(i, len(f1_gamete_list)):
            h_i = f1_gamete_list[i]
            h_j = f1_gamete_list[j]
            
            # The gametes produced by an F2 individual (h_i / h_j)
            recombinant_products = recombine(h_i, h_j)
            
            # Add these new gametes to the total set of F3 haplotypes
            f3_haplotypes.update(recombinant_products)

    # Step 4: Count and classify the unique sequences
    counts_by_transitions = collections.defaultdict(int)
    for h in sorted(list(f3_haplotypes)):
        t = count_transitions(h)
        counts_by_transitions[t] += 1
        
    print("The unique sequences found in the F3 generation can be classified by their number of allele transitions:")
    
    equation_parts = []
    for t in sorted(counts_by_transitions.keys()):
        count = counts_by_transitions[t]
        print(f"Number of sequences with {t} transition(s): {count}")
        equation_parts.append(str(count))
        
    total_sequences = len(f3_haplotypes)
    equation_str = " + ".join(equation_parts)
    
    print("\nThe total number of possible unique sequences is the sum of these counts:")
    print(f"Total = {equation_str} = {total_sequences}")

# Run the simulation and print the results
solve_genetics_problem()
print("\n<<<22>>>")