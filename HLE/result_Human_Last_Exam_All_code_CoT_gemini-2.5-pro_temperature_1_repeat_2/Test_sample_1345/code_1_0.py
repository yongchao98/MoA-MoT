import sys

def solve_for_n(N):
    """
    Calculates the maximal possible number of complex zeros for a given N.
    A zero k_1 is complex if the condition k_j is not real is met.
    """
    if not isinstance(N, int) or N <= 0:
        print("Error: N must be a positive integer.", file=sys.stderr)
        return

    # The case N=1 is special.
    # The governing polynomial has degree 1, so its single root must be real.
    # Thus, there are no complex zeros.
    if N == 1:
        result = 0
        factor = "N/A"
        base = "N/A"
        power = "N/A"
        print("For N = 1, the maximal number of complex zeros is 0.")

    # For N >= 2, the degree of the characteristic polynomial is N * 2^(N-1).
    # This degree is even, so all roots can be complex.
    else:
        factor = N
        base = 2
        power = N - 1
        result = factor * (base**power)
        # We print the components of the equation for clarity.
        print(f"For N = {N}, the maximal number of complex zeros is given by the formula:")
        print(f"Max Zeros = N * 2^(N-1)")
        print(f"Max Zeros = {factor} * {base}^({power})")
        print(f"Result = {result}")

def main():
    """
    Main function to parse input and call the solver.
    """
    # Use a default example value if no command-line argument is given.
    if len(sys.argv) > 1:
        try:
            N = int(sys.argv[1])
        except ValueError:
            print(f"Error: Invalid input '{sys.argv[1]}'. Please provide an integer.", file=sys.stderr)
            sys.exit(1)
    else:
        # Let's use N=3 as an example.
        print("No value for N provided. Using N=3 as an example.")
        N = 3

    solve_for_n(N)


if __name__ == "__main__":
    main()