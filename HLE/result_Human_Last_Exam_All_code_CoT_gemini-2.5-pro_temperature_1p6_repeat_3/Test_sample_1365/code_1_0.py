import argparse
import sys

def solve_experts_problem_bound(n, c):
    """
    Calculates and prints the upper bound on the number of mistakes made by
    a majority voting algorithm in a variant of the experts problem.

    The derivation for the upper bound M is:
    M <= (c * (n - 1) / 2) + (m_true)
    where m_true is the number of mistakes by the true expert.

    To get the upper bound, we use the maximum possible value for m_true,
    which is c - 1.
    M_bound = (c * (n - 1) / 2) + (c - 1)

    Args:
        n (int): The total number of experts.
        c (int): The mistake threshold for removing a false expert.
    """
    if n < 1 or c < 1:
        print("Error: Number of experts (n) and mistake threshold (c) must be at least 1.", file=sys.stderr)
        return

    # Using floating point division for generality
    term1 = (c * (n - 1)) / 2.0
    term2 = c - 1
    
    bound = term1 + term2

    print("--- Upper Bound Calculation ---")
    print(f"Given n = {n} experts and a mistake threshold c = {c}.")
    print("\nThe formula for the upper bound (M) on algorithm mistakes is:")
    print("M <= (c * (n - 1) / 2) + (c - 1)")
    
    print("\nSubstituting the given values:")
    # The final code should output each number in the final equation
    print(f"M <= ({c} * ({n} - 1) / 2) + ({c} - 1)")
    print(f"M <= ({c} * {n - 1} / 2) + {term2}")
    print(f"M <= ({c * (n - 1)} / 2) + {term2}")
    print(f"M <= {term1} + {term2}")
    
    # Check if the result is an integer for cleaner display
    if bound == int(bound):
        print(f"\nThe calculated upper bound is: {int(bound)}")
    else:
        print(f"\nThe calculated upper bound is: {bound}")

def main():
    """Main function to parse arguments and run the calculation."""
    parser = argparse.ArgumentParser(
        description="Calculate the upper bound on mistakes in a variant of the experts problem. Please provide values for n and c.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-n', 
        type=int, 
        required=True, 
        help="The total number of experts (e.g., 10)."
    )
    parser.add_argument(
        '-c', 
        type=int, 
        required=True, 
        help="The mistake threshold for removing an expert (e.g., 5)."
    )
    
    # If no arguments are provided, print help and exit.
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    solve_experts_problem_bound(args.n, args.c)

if __name__ == "__main__":
    main()