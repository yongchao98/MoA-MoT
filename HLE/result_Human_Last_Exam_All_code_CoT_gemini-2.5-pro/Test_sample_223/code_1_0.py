def solve_chair_puzzle():
    """
    Solves the chair puzzle by deriving and applying a mathematical inequality.
    """
    # Total number of chairs
    n = 20

    # Let k be the maximum number of occupied chairs.

    # Step 1: Explain the logic for maximizing occupants.
    # The number of occupants increases only when a person sits in a chair with no neighbors.
    # To maximize the final count, the strategy must be to always choose such a chair.
    # This implies that no two occupied chairs can be adjacent.

    # Step 2: Formulate the mathematical relationship.
    # If k chairs are occupied by non-adjacent people, they must be separated by at least one empty chair.
    # A minimal arrangement would be X_X_X_..._X.
    # This arrangement has k occupied chairs ('X') and k-1 empty chairs ('_') separating them.
    # The minimum number of chairs required for such an arrangement is k + (k-1), which is 2*k - 1.

    # Step 3: Create and solve the inequality.
    # The minimum chairs required (2*k - 1) must be less than or equal to the total chairs available (n).
    # This gives us the inequality: 2*k - 1 <= n

    print("To find the maximum number of occupied chairs, we follow this logic:")
    print("1. The number of people can only increase if a new person sits in a chair with no occupied neighbors.")
    print("2. Therefore, to maximize the count, all occupied chairs must be non-adjacent.")
    print("\nLet k be the number of occupied chairs and n be the total number of chairs.")
    print("For k people to be non-adjacent, we need at least k chairs for them and k-1 empty chairs as separators.")
    print("This leads to the inequality: 2*k - 1 <= n")

    # Step 4: Substitute the value of n and calculate k.
    print("\nNow, we solve for k with n = 20:")
    
    # The following print statements explicitly show each number in the equation.
    print(f"2 * k - 1 <= {n}")
    
    n_plus_1 = n + 1
    print(f"2 * k <= {n} + 1")
    print(f"2 * k <= {n_plus_1}")
    
    k_max_float = n_plus_1 / 2
    print(f"k <= {n_plus_1} / 2")
    print(f"k <= {k_max_float}")

    # Since k must be an integer, we take the floor of the result.
    k_max_int = int(k_max_float)

    print("\nSince the number of people must be an integer, the maximum is the whole number part of the result.")
    print(f"The maximum number of chairs that can be occupied is {k_max_int}.")

solve_chair_puzzle()