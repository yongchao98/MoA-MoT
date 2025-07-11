def solve_tree_diameter():
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.
    """
    try:
        input_str = input("Enter the positive integers n and m, separated by a space: ")
        n_str, m_str = input_str.split()
        n = int(n_str)
        m = int(m_str)

        # A tree must have n>=1 (for n+2>=3 vertices) and at least 2 leaves.
        # The maximum number of leaves in a tree with n+2 vertices is n+1 (a star graph).
        if not (n >= 1 and 2 <= m <= n + 1):
             print(f"\nInvalid input combination.")
             print(f"Constraints: n must be a positive integer, and for n={n}, m must be between 2 and {n+1}, inclusive.")
             return

        # I is the number of internal nodes (degree > 1)
        i_val = n + 2 - m

        # If there are many internal nodes (I>=4) and enough leaves to make them
        # form a compact star-like structure (m >= I-1), the diameter is a small constant.
        if i_val >= 4 and m >= i_val - 1:
            diameter = 4
            print("\nThe number of internal nodes I = n + 2 - m = {} + 2 - {} = {}.".format(n, m, i_val))
            print("The conditions for a compact 'star-of-internals' structure are met (I >= 4 and m >= I-1).")
            print("The minimum possible diameter is 4.")
        else:
            # Otherwise, the minimum diameter is achieved with a 'caterpillar' structure,
            # where the diameter is determined by the length of the internal node path.
            diameter = n + 3 - m
            print("\nThe minimum diameter is given by the formula d = n + 3 - m.")
            # Printing each number in the final equation, as requested.
            print(f"{n} + 3 - {m} = {diameter}")

    except (ValueError, IndexError):
        print("\nInvalid input. Please enter two integers separated by a space.")

solve_tree_diameter()