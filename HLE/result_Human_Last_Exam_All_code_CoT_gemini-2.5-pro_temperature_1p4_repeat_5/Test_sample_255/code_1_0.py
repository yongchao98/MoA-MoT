import math

def phi(n):
    """
    Computes Euler's totient function phi(n).
    For n = 2^k, phi(n) = n - n/2 = n/2.
    """
    if n == 1:
        return 1
    # A simple check for powers of 2, which is all we need here.
    if n > 0 and (n & (n-1) == 0):
        return n // 2
    
    # General case implementation for completeness.
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def solve_cohomology_dimension():
    """
    Calculates the dimension of H^2(G,M) based on the derived formula.
    The dimension is the sum of phi(d) for d in {2, 4, 8}.
    """
    
    print("Based on the analysis, the dimension of H^2(G,M) is the sum of contributions from submodules M_d for d dividing 8.")

    # d values for which the dimension contribution is non-zero
    d_values_nonzero = [2, 4, 8]
    
    # Calculate dimensions for each d that contributes
    dimensions_to_sum = []
    
    # Add the dimension for d=1, which is 0
    dim_d1 = 0
    print(f"Contribution from d=1: {dim_d1}")
    dimensions_to_sum.append(dim_d1)
    
    # Calculate dimensions for d=2, 4, 8
    for d in d_values_nonzero:
        dim_d = phi(d)
        print(f"Contribution from d={d}: phi({d}) = {dim_d}")
        dimensions_to_sum.append(dim_d)
        
    total_dimension = sum(dimensions_to_sum)
    
    print("\nThe other d values (16, 32, 64, 128) that divide 128 have a zero contribution.")
    
    # Build and print the final equation string as requested.
    equation_parts = [str(dim) for dim in dimensions_to_sum]
    final_equation = " + ".join(equation_parts) + f" = {total_dimension}"
    
    print("\nThe total dimension is the sum of these contributions:")
    print(final_equation)

solve_cohomology_dimension()