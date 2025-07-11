def solve_coin_puzzle():
    """
    Calculates the maximum number of real coins that can be guaranteed
    to be identified in two weighings.
    """

    # Let n be the number of coins on each pan in the first weighing.
    # There are two main scenarios for the first weighing's outcome.

    # Scenario 1: The scale is unbalanced.
    # The guaranteed number of real coins is 1000 - n.
    # Let's represent this with a function.
    def unbalanced_guarantee(n):
        return 1000 - n

    # Scenario 2: The scale is balanced.
    # An optimal second weighing can prove that the 2n coins from the
    # first weighing are real.
    # The guaranteed number of real coins is 2n.
    def balanced_guarantee(n):
        return 2 * n

    # To find the maximum number of coins we can *guarantee*, we need to
    # maximize the minimum of these two outcomes.
    # We find the approximate n where the outcomes are equal:
    # 1000 - n = 2n  => 3n = 1000 => n = 333.33...
    
    # Since n must be an integer, we check the integers around this value.
    n1 = 333
    n2 = 334
    
    # Calculate the guaranteed number of coins for n = 333
    guarantee_n1 = min(unbalanced_guarantee(n1), balanced_guarantee(n1))
    
    # Calculate the guaranteed number of coins for n = 334
    guarantee_n2 = min(unbalanced_guarantee(n2), balanced_guarantee(n2))

    # The maximum guaranteed number is the maximum of these results.
    max_guaranteed_coins = max(guarantee_n1, guarantee_n2)
    
    # We print the logic for n=334 as an example.
    n = 334
    print(f"Let n = {n} be the number of coins on each pan in the first weighing.")
    print("If the scale is unbalanced, the number of guaranteed real coins is:")
    print(f"1000 - n = 1000 - {n} = {1000-n}")
    print("If the scale is balanced, an optimal second weighing can guarantee:")
    print(f"2 * n = 2 * {n} = {2*n}")
    print("The number of coins we can guarantee to identify is the minimum of these two outcomes.")
    print(f"min({1000-n}, {2*n}) = {max_guaranteed_coins}")
    print("\nFinal Answer:")
    print(max_guaranteed_coins)


solve_coin_puzzle()