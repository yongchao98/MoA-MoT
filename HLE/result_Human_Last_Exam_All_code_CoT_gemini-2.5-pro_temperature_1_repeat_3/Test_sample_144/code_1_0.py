import collections

def solve_genetics_sequences():
    """
    Calculates the number of unique autosome sequences in the F3 generation
    based on a specific model of recombination.
    """
    num_snps = 5

    # Step 1: Define parental chromosomes
    p1_chr = 'A' * num_snps
    p2_chr = 'B' * num_snps
    
    # Step 2: Generate F1 gametes from the F1 individual (p1_chr / p2_chr)
    # A single crossover occurs at one of the (num_snps - 1) intervals.
    s_f1_gametes = set()
    for k in range(1, num_snps):
        s_f1_gametes.add(p1_chr[:k] + p2_chr[k:])
        s_f1_gametes.add(p2_chr[:k] + p1_chr[k:])

    # Step 3: Generate F2 gametes (these form the F3 generation)
    # F2 individuals are formed by pairing any two F1 gametes.
    s_f2_gametes = set()
    f1_gametes_list = list(s_f1_gametes)

    for c1 in f1_gametes_list:
        for c2 in f1_gametes_list:
            # An F2 individual has genotype c1 / c2.
            # Generate gametes from this individual via single crossover.
            for k in range(1, num_snps):
                s_f2_gametes.add(c1[:k] + c2[k:])
                s_f2_gametes.add(c2[:k] + c1[k:])

    # Step 4: Classify the unique F2 gametes by number of transitions
    def count_transitions(seq):
        """Counts the number of allele changes in a sequence."""
        transitions = 0
        for i in range(len(seq) - 1):
            if seq[i] != seq[i+1]:
                transitions += 1
        return transitions

    transition_counts = collections.defaultdict(int)
    for seq in s_f2_gametes:
        t = count_transitions(seq)
        transition_counts[t] += 1
        
    # Step 5: Print the results, including the "equation"
    print("The unique sequences found in the F3 generation can be classified by their number of allele transitions:")
    
    equation_parts = []
    for t in sorted(transition_counts.keys()):
        count = transition_counts[t]
        print(f"- Sequences with {t} transitions: {count}")
        equation_parts.append(str(count))
        
    final_count = len(s_f2_gametes)
    equation_str = " + ".join(equation_parts)
    
    print("\nThe final count is the sum of these numbers:")
    print(f"{equation_str} = {final_count}")

solve_genetics_sequences()
<<<26>>>