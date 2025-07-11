def solve_coin_puzzle():
    """
    Calculates the maximum number of real coins that can be guaranteed to be
    identified in two weighings.

    The logic is based on balancing the number of identified coins from each
    possible outcome of the weighings.

    Let n be the number of coins on each side of the scale in the first weighing.
    Let m be the number of coins on each side of the scale in the second weighing.

    The number of guaranteed real coins in each of the three main scenarios is:
    1. First weighing is imbalanced: 1000 - n
    2. First weighing is balanced, second is imbalanced: 1000 - m
    3. First and second weighings are balanced: 2*n + 2*m

    To maximize the minimum guaranteed number, we set these equal to each other:
    1000 - n = 1000 - m  => n = m
    1000 - n = 2*n + 2*m

    Substituting n = m into the second equation:
    1000 - n = 2*n + 2*n
    1000 - n = 4*n
    1000 = 5*n
    n = 1000 / 5
    """

    # Total coins
    total_coins = 1000
    
    # We solve the equation 1000 = 5 * n to find n
    # The equation can be written as:
    # 5 * n = 1000
    n = total_coins / 5
    m = n
    
    guaranteed_real_coins = total_coins - n
    
    print("The strategy is as follows:")
    print(f"1. First Weighing: Place {int(n)} coins on the left and {int(n)} coins on the right.")
    print(f"   - If they are imbalanced, the {int(1000-n)} coins not on the lighter side are real.")
    print(f"   - If they are balanced, the {int(2*n)} coins on the scale are real. The 4 fakes are in the remaining {int(1000-2*n)} coins.")
    print(f"2. Second Weighing (only if first was balanced): From the {int(1000-2*n)} remaining coins, weigh {int(m)} vs {int(m)}.")
    print(f"   - If imbalanced, total real coins = {int(2*n)} (from W1) + {int(m)} + {int(1000-2*n-2*m)} (from W2) = {int(1000-m)}.")
    print(f"   - If balanced, total real coins = {int(2*n)} (from W1) + {int(2*m)} (from W2) = {int(2*n + 2*m)}.")
    
    print("\nTo guarantee the maximum number of real coins, we need the outcome to be the same in all cases.")
    print(f"Equation: 1000 - n = 2 * n + 2 * n")
    print(f"1000 = 5 * n")
    print(f"n = {int(n)}")
    
    print("\nThis strategy guarantees identifying the same number of real coins regardless of the outcome:")
    print(f"Case 1 result: 1000 - {int(n)} = {int(guaranteed_real_coins)}")
    print(f"Case 2a result: 1000 - {int(m)} = {int(guaranteed_real_coins)}")
    print(f"Case 2b result: 2 * {int(n)} + 2 * {int(m)} = {int(guaranteed_real_coins)}")
    
    print("\nMaximum number of real coins you can guarantee to identify:")
    print(int(guaranteed_real_coins))

solve_coin_puzzle()