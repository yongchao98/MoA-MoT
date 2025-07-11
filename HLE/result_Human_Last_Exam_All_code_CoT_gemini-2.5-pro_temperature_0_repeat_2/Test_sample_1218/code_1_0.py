def solve_max_n(k):
    """
    Calculates the maximum value of n in terms of k for a k-uniform
    intersecting family with full differences of size k-1.

    The formula derived from combinatorial arguments is n = k^2 - k + 1.

    Args:
        k (int): The size of the subsets in the family (must be >= 2).
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.")
        return

    # Calculate the maximum value of n using the formula n = k^2 - k + 1
    n_max = k**2 - k + 1

    # Output the final equation and the result
    print(f"For k = {k}, the maximum value of n is given by the equation:")
    print(f"n = k^2 - k + 1")
    print(f"n = {k}^2 - {k} + 1")
    print(f"n = {k*k} - {k} + 1")
    print(f"n = {n_max}")

if __name__ == '__main__':
    # Example usage:
    # Get k from the user
    try:
        k_input = int(input("Enter the value of k (an integer >= 2): "))
        solve_max_n(k_input)
    except ValueError:
        print("Invalid input. Please enter an integer.")
