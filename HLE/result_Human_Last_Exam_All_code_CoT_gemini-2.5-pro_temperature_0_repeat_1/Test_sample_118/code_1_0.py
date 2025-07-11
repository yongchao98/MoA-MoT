def solve_coin_puzzle():
    """
    This function calculates the maximum number of real coins that can be guaranteed to be identified.
    
    Problem statement:
    - 1000 total coins
    - 4 fake coins (lighter)
    - 996 real coins
    - 2 weighings on a balance scale

    The solution to this problem is based on isolating a small group of coins
    that contains all the fakes, thereby guaranteeing the remaining coins are real.
    It can be shown that with 2 weighings, we can isolate the 4 fakes into a group of 8 coins.
    This means the remaining 1000 - 8 = 992 coins can be guaranteed as real.
    """
    
    total_coins = 1000
    
    # The number of coins that we cannot guarantee to be real consists of the 4 fake coins
    # plus 4 real coins that get mixed in with the fakes in the ambiguous group.
    # This group of 8 coins contains all 4 fakes.
    unidentified_group_size = 8
    
    guaranteed_real_coins = total_coins - unidentified_group_size
    
    print("Total coins: 1000")
    print("Fake coins: 4 (lighter)")
    print("Weighings allowed: 2")
    print("The strategy involves isolating the 4 fake coins into a small group of 8 coins.")
    print("This guarantees the remaining coins are all real.")
    print(f"The number of guaranteed real coins is Total Coins - Size of Unidentified Group.")
    print(f"{total_coins} - {unidentified_group_size} = {guaranteed_real_coins}")

solve_coin_puzzle()