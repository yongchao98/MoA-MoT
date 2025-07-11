import sys

def solve_max_zeros():
    """
    Calculates the maximal number of complex zeros for the given matrix determinant equation.

    The problem asks for the maximal number of complex zeros for det(B(k))=0, where
    B(k) is an N x N matrix. The number of such zeros is given by the formula N * 2^(N-1).
    This script calculates and prints this value for a given N.
    """
    try:
        # You can change the value of N here
        N = 3
        # For command line execution, you can provide N as an argument
        if len(sys.argv) > 1:
            N = int(sys.argv[1])
        
        if not isinstance(N, int) or N <= 0:
            raise ValueError("N must be a positive integer.")

        power = N - 1
        base = 2
        
        # Using integer exponentiation
        result = N * (base ** power)
        
        print(f"For N = {N}:")
        print(f"The formula for the maximal number of complex zeros is N * 2^(N-1).")
        # Outputting each number in the final equation
        print(f"The calculation is: {N} * {base}^({N} - 1) = {result}")

    except (ValueError, IndexError) as e:
        print(f"Error: {e}", file=sys.stderr)
        print("Please provide a single positive integer for N.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    solve_max_zeros()
