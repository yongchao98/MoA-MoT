def solve_autosome_sequences():
    """
    Calculates the number of unique autosome sequences in the F3 generation
    based on the specified genetic crossing and recombination rules.
    """

    # The number of SNPs on the autosome.
    num_snps = 5
    parent1 = "A" * num_snps
    parent2 = "B" * num_snps
    
    # === F1 Gametes (Haplotypes in F2 Generation) ===
    # With "exactly one" crossover per gamete, only recombinant haplotypes are formed.
    # There are (num_snps - 1) locations for a crossover. Each creates 2 unique sequences.
    f1_gametes = set()
    num_crossover_points = num_snps - 1

    for k in range(1, num_crossover_points + 1):
        f1_gametes.add(parent1[:k] + parent2[k:])
        f1_gametes.add(parent2[:k] + parent1[k:])

    # === F2 Gametes (Haplotypes in F3 Generation) ===
    # An F2 individual's genotype is (h1/h2) where h1 and h2 are from f1_gametes.
    # We find all possible gametes from all possible F2 genotypes.
    f2_gametes = set()
    f2_possible_haplotypes = list(f1_gametes)

    for h1 in f2_possible_haplotypes:
        for h2 in f2_possible_haplotypes:
            for k in range(1, num_crossover_points + 1):
                # Generate the two reciprocal recombinant gametes
                f2_gametes.add(h1[:k] + h2[k:])
                f2_gametes.add(h2[:k] + h1[k:])
    
    # === Categorize and Count Final Sequences ===
    # We can categorize the final sequences by the number of "switches" (recombination events in their history).
    seqs_by_switch_count = {i: 0 for i in range(num_snps)}
    for seq in f2_gametes:
        switches = 0
        for i in range(num_snps - 1):
            if seq[i] != seq[i+1]:
                switches += 1
        if switches < num_snps:
            seqs_by_switch_count[switches] += 1
            
    num_0_switch = seqs_by_switch_count[0]
    num_1_switch = seqs_by_switch_count[1]
    num_2_switch = seqs_by_switch_count[2]
    total_sequences = len(f2_gametes)

    print("The number of possible unique sequences can be calculated by summing the types:")
    print(f"Parental-like sequences (0 switches): {num_0_switch}")
    print(f"Sequences with 1 switch point: {num_1_switch}")
    print(f"Sequences with 2 switch points: {num_2_switch}")
    print("\nFinal calculation:")
    print(f"{num_0_switch} + {num_1_switch} + {num_2_switch} = {total_sequences}")

solve_autosome_sequences()
<<<22>>>