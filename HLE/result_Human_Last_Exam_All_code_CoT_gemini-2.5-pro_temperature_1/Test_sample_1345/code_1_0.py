import sys

def calculate_max_complex_zeros(N):
    """
    Calculates the maximal possible number of complex zeros for a given matrix size N.

    The maximal number of complex zeros is given by the degree of a polynomial
    derived from the determinant equation. The formula for this degree (for N >= 1) is:
    D_N = 2^(N+1) - 4

    Args:
        N (int): The dimension of the matrix, must be a positive integer.

    Returns:
        None. The function prints the result directly.
    """
    if not isinstance(N, int) or N < 1:
        print("Error: N must be a positive integer.")
        return

    # The formula is valid for N >= 1.
    power = N + 1
    # In Python, ** has precedence over -, so parentheses are not strictly needed
    # but are kept for clarity.
    result = (2**power) - 4
    
    print(f"For a matrix of size N = {N}:")
    print("The maximal possible number of complex zeros is given by the formula 2^(N+1) - 4.")
    print("\nThe calculation is as follows:")
    # Showing each number in the final equation as requested
    print(f"2^({N} + 1) - 4 = 2^{power} - 4 = {2**power} - 4 = {result}")

def main():
    """
    Main function to parse the command-line argument for N and run the calculation.
    """
    if len(sys.argv) != 2:
        print("Usage: python", sys.argv[0], "N")
        print("  N: The dimension of the N x N matrix (a positive integer).")
        sys.exit(1)

    try:
        N_arg = int(sys.argv[1])
        calculate_max_complex_zeros(N_arg)
    except ValueError:
        print("Error: The provided argument for N must be an integer.")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()