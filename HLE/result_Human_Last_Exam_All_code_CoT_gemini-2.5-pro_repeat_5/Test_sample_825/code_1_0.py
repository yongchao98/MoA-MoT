def count_distinct_polynomials():
    """
    This script calculates the number of distinct polynomials p(n) that can occur
    as the dimension of an FS_n-submodule of V_n.

    The method is to count all possible combinations of submodule dimensions
    and then subtract the redundancies caused by linear dependencies among the
    basis dimension polynomials.
    """

    # The dimension p(n) of a submodule is a sum:
    # p(n) = c1*dim(M1) + c2*dim(M2) + c3*dim(M3) + c4*dim(M4)
    # The number of choices for the coefficients c_i depends on the multiplicities
    # of the irreducible components M_i in V_n for n >= 4.
    # V_n decomposes as 2*M1 + 2*M2 + M3 + M4.

    # Number of choices for c1 from the 2*M1 component (submodules 0, M1, 2*M1)
    c1_choices = 3
    # Number of choices for c2 from the 2*M2 component (submodules 0, M2, 2*M2)
    c2_choices = 3
    # Number of choices for c3 from the M3 component (submodules 0, M3)
    c3_choices = 2
    # Number of choices for c4 from the M4 component (submodules 0, M4)
    c4_choices = 2

    # Total number of combinations of coefficients
    total_combinations = c1_choices * c2_choices * c3_choices * c4_choices

    print("Step 1: Calculate the total number of coefficient combinations.")
    print(f"Number of choices for c1, c2, c3, c4 are {c1_choices}, {c2_choices}, {c3_choices}, {c4_choices} respectively.")
    print(f"Total combinations = {c1_choices} * {c2_choices} * {c3_choices} * {c4_choices} = {total_combinations}")
    print("-" * 20)

    # A linear dependency dim(M1) + dim(M3) = dim(M2) + dim(M4) causes redundancy.
    # We find the number of pairs of coefficient vectors that give the same polynomial.
    # A pair exists if C = (c1, c2, c3, c4) and C' = (c1-1, c2+1, c3-1, c4+1) are both valid.
    # This happens if c1>=1, c2<=1, c3=1, c4=0.

    # Number of choices for c1 in {1, 2}
    source_c1_choices = 2
    # Number of choices for c2 in {0, 1}
    source_c2_choices = 2
    # Number of choices for c3 must be 1
    source_c3_choices = 1
    # Number of choices for c4 must be 0
    source_c4_choices = 1

    redundant_pairs = source_c1_choices * source_c2_choices * source_c3_choices * source_c4_choices

    print("Step 2: Calculate the number of redundant pairs due to the dependency.")
    print(f"Number of redundant pairs = {source_c1_choices} * {source_c2_choices} * {source_c3_choices} * {source_c4_choices} = {redundant_pairs}")
    print("-" * 20)

    # The final number of distinct polynomials is the total minus the redundant pairs.
    distinct_polynomials = total_combinations - redundant_pairs

    print("Step 3: Calculate the final number of distinct polynomials.")
    print("The final equation is: Distinct Polynomials = Total Combinations - Redundant Pairs")
    print(f"The result is: {distinct_polynomials} = {total_combinations} - {redundant_pairs}")

# Run the calculation
count_distinct_polynomials()