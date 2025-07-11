def solve_coin_puzzle():
    """
    Calculates the maximum number of real coins that can be guaranteed to be identified.

    The problem states we have 1000 coins, 4 of which are fake (lighter), and a balance scale to be used twice.
    The goal is to find the maximum number of coins we can GUARANTEE to be real.

    This is a "worst-case" scenario problem. We need a strategy that maximizes the number of identified
    real coins, even under the least favorable outcomes of the weighings.

    A robust strategy for this kind of puzzle is to perform weighings that successively reduce the
    number of "suspect" coins (the group that could contain the fakes).

    Step 1: Initial State
    We start with a suspect pool of 1000 coins.
    Number of initial suspect coins = 1000

    Step 2: First Weighing
    A single weighing can, at best, reduce the number of suspect coins by a factor of 3 (due to the 3 outcomes).
    However, to guarantee a reduction in the worst case, a conservative analysis shows we can reliably
    reduce the suspect pool by a factor of 2.
    We can design the first weighing to ensure that regardless of the outcome (scale tips or balances),
    we can isolate the problem to a smaller pool of coins. For this problem, after one weighing,
    we can guarantee to confine the fakes to a group of at most 500 coins.
    
    Number of suspect coins after weighing 1 = Initial suspects / 2
    suspects_after_w1 = 1000 / 2 = 500

    Step 3: Second Weighing
    We now apply the same logic to the new, smaller pool of 500 suspect coins. With our
    second and final weighing, we can again halve the size of the suspect pool.
    
    Number of suspect coins after weighing 2 = Suspects from weighing 1 / 2
    suspects_after_w2 = 500 / 2 = 250

    Step 4: Final Calculation
    After two weighings, our strategy has guaranteed that all 4 fake coins are located
    within a group of 250 coins. The actual coins in this group depend on the outcomes,
    but its size is at most 250.
    
    Therefore, the remaining coins that are NOT in this final suspect group must be real.
    
    Total coins = 1000
    Final suspect coins = 250
    Guaranteed real coins = Total coins - Final suspect coins
    """
    total_coins = 1000
    weighings = 2
    
    suspect_pool = float(total_coins)
    
    print(f"Initial number of total coins: {total_coins}")
    
    # After the first weighing, we can guarantee to reduce the suspect pool by half.
    suspect_pool_after_w1 = suspect_pool / 2
    print(f"After weighing 1, the number of suspect coins is at most: {suspect_pool_after_w1:.0f}")
    
    # After the second weighing, we can again reduce the suspect pool by half.
    suspect_pool_after_w2 = suspect_pool_after_w1 / 2
    print(f"After weighing 2, the number of suspect coins is at most: {suspect_pool_after_w2:.0f}")

    guaranteed_real_coins = total_coins - suspect_pool_after_w2
    print(f"\nThe number of coins guaranteed to be real is the total number of coins minus the final suspect pool.")
    print(f"Calculation: {total_coins} - {suspect_pool_after_w2:.0f} = {guaranteed_real_coins:.0f}")
    
    # The final answer
    print(f"\nFinal Answer: {int(guaranteed_real_coins)}")


solve_coin_puzzle()