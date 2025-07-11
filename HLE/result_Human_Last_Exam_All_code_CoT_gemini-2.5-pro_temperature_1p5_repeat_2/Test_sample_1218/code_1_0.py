import sys

def solve_combinatorial_problem(k):
    """
    Calculates the maximum value of n in terms of k for a k-uniform intersecting family
    with full differences of size k-1.

    The formula is derived from the properties of projective planes, which provide
    a valid construction for k=2 and k=3. The proposed maximum value is n = k^2 - k + 1.

    Args:
        k (int): The size of the subsets in the family. Must be an integer >= 2.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.", file=sys.stderr)
        return

    # The maximum value of n is given by the formula k^2 - k + 1.
    n = k**2 - k + 1

    # Output the result in a clear format.
    print(f"For k = {k}, the maximum value of n is calculated as:")
    print(f"n = k^2 - k + 1")
    print(f"n = {k}^2 - {k} + 1")
    print(f"n = {k*k} - {k} + 1")
    print(f"n = {n}")


if __name__ == '__main__':
    # Example usage of the function.
    # You can change the value of k here.
    try:
        if len(sys.argv) > 1:
            k_value = int(sys.argv[1])
        else:
            # Default k value if not provided as a command-line argument.
            print("Usage: python your_script_name.py <k_value>")
            print("Using default k=4 for demonstration.")
            k_value = 4
        solve_combinatorial_problem(k_value)
    except (ValueError, IndexError):
        print("Error: Please provide a valid integer for k.", file=sys.stderr)
