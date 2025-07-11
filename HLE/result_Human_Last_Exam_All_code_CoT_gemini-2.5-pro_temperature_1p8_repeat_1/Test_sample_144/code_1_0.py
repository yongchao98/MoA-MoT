import collections

def solve_genetics_recombination():
    """
    Calculates the number of unique autosome sequences in the F3 generation.
    """
    num_snps = 5

    # Step 1: Define parental chromosomes
    parent_a = 'A' * num_snps
    parent_b = 'B' * num_snps

    def get_all_single_recombinants(chrom1, chrom2):
        """Generates a set of all single-crossover recombinants."""
        length = len(chrom1)
        recombinants = set()
        # There are (length - 1) possible crossover points
        for i in range(1, length):
            recombinant1 = chrom1[:i] + chrom2[i:]
            recombinants.add(recombinant1)
            recombinant2 = chrom2[:i] + chrom1[i:]
            recombinants.add(recombinant2)
        return recombinants

    # Step 2: Generate the 8 unique gametes from the F1 generation
    f1_gametes = get_all_single_recombinants(parent_a, parent_b)

    # Step 3: Generate all possible gametes from the F2 generation.
    # An F2 individual's genotype is a pair of chromosomes from f1_gametes.
    f2_gametes = set()
    f1_gametes_list = list(f1_gametes)
    for i in range(len(f1_gametes_list)):
        for j in range(i, len(f1_gametes_list)):
            c1 = f1_gametes_list[i]
            c2 = f1_gametes_list[j]
            # Generate gametes from this F2 individual and add to the set
            recombinants_from_f2 = get_all_single_recombinants(c1, c2)
            f2_gametes.update(recombinants_from_f2)

    # Step 4: Categorize the resulting unique sequences by number of switches
    def count_switches(seq):
        """Counts the number of allele type switches in a sequence."""
        switches = 0
        for i in range(len(seq) - 1):
            if seq[i] != seq[i+1]:
                switches += 1
        return switches

    switch_counts = collections.defaultdict(int)
    for seq in f2_gametes:
        s = count_switches(seq)
        switch_counts[s] += 1

    num_0_switch = switch_counts.get(0, 0)
    num_1_switch = switch_counts.get(1, 0)
    num_2_switch = switch_counts.get(2, 0)
    
    total = num_0_switch + num_1_switch + num_2_switch

    print("The unique sequences in the F3 generation can be categorized by the number of switches between 'A' and 'B' alleles.")
    print(f"Number of sequences with 0 switches (e.g., AAAAA): {num_0_switch}")
    print(f"Number of sequences with 1 switch (e.g., AABBB): {num_1_switch}")
    print(f"Number of sequences with 2 switches (e.g., AABAA): {num_2_switch}")
    print("No sequences with more than 2 switches are formed in this generation.")
    print("\nThe total number of possible unique sequences is the sum of these counts:")
    print(f"{num_0_switch} + {num_1_switch} + {num_2_switch} = {total}")

solve_genetics_recombination()