import sys

def solve_max_n(k):
    """
    This function calculates the maximum value of n for a given k based on the
    combinatorial properties described.

    The derivation shows that n must be at least 2k - 1.
    A construction exists for n = 2k - 1, proving it is also the maximum.

    Args:
        k (int): The size of the uniform subsets in the intersecting family.
    """
    if not isinstance(k, int) or k < 1:
        print("Error: k must be a positive integer.", file=sys.stderr)
        return

    # Based on the derivation, the maximum value of n is 2k - 1.
    n = 2 * k - 1

    print(f"For k = {k}, the maximum value of n is determined by the equation:")
    print(f"n_max = 2 * k - 1")
    print(f"n_max = 2 * {k} - 1")
    print(f"n_max = {n}")


if __name__ == '__main__':
    # You can run this script with a command-line argument for k.
    # Example: python your_script_name.py 5
    if len(sys.argv) > 1:
        try:
            k_from_arg = int(sys.argv[1])
            solve_max_n(k_from_arg)
        except ValueError:
            print("Invalid input. Please provide an integer for k.", file=sys.stderr)
    else:
        # If no argument is provided, run with an example value.
        print("No command-line argument for k provided. Using k=4 as an example.")
        k_example = 4
        solve_max_n(k_example)
