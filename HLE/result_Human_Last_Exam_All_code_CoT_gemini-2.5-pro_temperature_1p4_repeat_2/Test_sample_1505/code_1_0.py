import sys

def calculate_sum_approximation(n):
    """
    Calculates the approximation for the sum for a given positive integer n.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.", file=sys.stderr)
        return

    # The approximation formula is: n^2/2 + 1/120 + 1/(252*n)
    term1 = n**2 / 2
    term2 = 1 / 120
    term3 = 1 / (252 * n)

    # The sum is the approximation
    result = term1 + term2 + term3

    # Print the terms of the equation and the final result
    print(f"The approximation is based on the formula: n^2/2 + 1/120 + 1/(252*n)")
    print(f"For n = {n}:")
    print(f"The term n^2/2 = {term1}")
    print(f"The term 1/120 = {term2}")
    print(f"The term 1/(252*n) = {term3}")
    print(f"The final approximated sum is: {result}")
    
    # Per instructions, also printing the final equation with values
    print("\nFinal equation with calculated values:")
    print(f"{term1} + {term2} + {term3} = {result}")

if __name__ == '__main__':
    # You can change this value to any positive integer.
    n_value = 100
    calculate_sum_approximation(n_value)
