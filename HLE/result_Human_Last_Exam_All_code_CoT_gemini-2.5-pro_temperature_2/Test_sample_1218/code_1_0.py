def solve_max_n(k):
    """
    Calculates the maximum value of n in terms of k based on the derivation.
    The formula is n = 2k - 1.

    Args:
        k (int): The size of the subsets in the family F. Must be an integer > 1.
    """
    if not isinstance(k, int) or k <= 1:
        print("Error: k must be an integer greater than 1.")
        return

    # Calculate n using the derived formula
    n = 2 * k - 1

    # Print the final equation with all its numbers
    print(f"For k = {k}, the maximum value of n is given by the equation:")
    # The final equation is n = 2*k - 1
    # We output each number involved in this equation.
    print(f"{n} = 2 * {k} - 1")

if __name__ == '__main__':
    # You can change the value of k here to test for different cases.
    # Let's use an example value for k, e.g., k = 4.
    example_k = 4
    solve_max_n(example_k)
    print("-" * 20)
    # Another example for k = 10
    example_k_2 = 10
    solve_max_n(example_k_2)
