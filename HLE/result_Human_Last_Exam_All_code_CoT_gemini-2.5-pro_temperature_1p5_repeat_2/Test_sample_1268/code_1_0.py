import sys

def calculate_norm_bound(N):
    """
    Calculates the upper bound for the maximum norm in relation to the covolume
    for a given squarefree natural number N.

    Args:
        N (int): A squarefree natural number, representing the discriminant.
    """
    if N <= 0:
        print("Error: The squarefree number N must be a positive integer.", file=sys.stderr)
        return

    # In the context of definite quaternion algebras, the covolume (V) of a
    # maximal order is related to the discriminant (N) by the formula: V = N / 2.
    covolume_V = N / 2.0

    # A result from the geometry of numbers states that there exists a basis for the
    # order where the norm of each basis element is bounded by N / 2.
    # So, the upper bound for the maximum norm (k) is k <= N / 2.
    # This means the upper bound is equal to the covolume V.
    upper_bound = covolume_V

    # Output the result, including the numbers in the final equation as requested.
    # We use 'k' to represent the maximum norm k_{k,inf}.
    print(f"For the squarefree number N = {N}:")
    print(f"The covolume V is calculated as N / 2 = {N} / 2 = {covolume_V}.")
    print("The relation for the upper bound on the maximum norm (k) is:")
    print(f"k <= {upper_bound}")

if __name__ == '__main__':
    # You can provide N as a command-line argument.
    # If no argument is given, a default value is used for demonstration.
    if len(sys.argv) > 1:
        try:
            input_N = int(sys.argv[1])
            calculate_norm_bound(input_N)
        except ValueError:
            print("Error: Please provide a valid integer for N.", file=sys.stderr)
            sys.exit(1)
    else:
        # Using a default squarefree number for demonstration, e.g., 3*5*7 = 105.
        print("No command-line argument provided for N. Using a default value N=105.")
        default_N = 105
        calculate_norm_bound(default_N)
