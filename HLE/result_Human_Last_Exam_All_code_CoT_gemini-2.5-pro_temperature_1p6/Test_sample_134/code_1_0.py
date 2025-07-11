import math

def solve():
    """
    Solves the coin game problem by analyzing the game's structure.
    """
    one_euro_coins = 136
    two_euro_coins = 87
    total_coins = one_euro_coins + two_euro_coins

    print(f"There are {one_euro_coins} 1-euro coins and {two_euro_coins} 2-euro coins.")
    print(f"The total number of coins in the line is {total_coins}.")
    print("")

    # Since the total number of coins is odd, the number of turns is not equal.
    player1_turns = math.ceil(total_coins / 2)
    player2_turns = math.floor(total_coins / 2)

    print("The game proceeds in alternating turns.")
    print(f"Player 1 gets to pick {player1_turns} coins in total.")
    print(f"Player 2 gets to pick {player2_turns} coins in total.")
    print("")
    
    print("Conclusion:")
    print("Player 1 gets one more coin than Player 2 ({0} vs {1}).".format(player1_turns, player2_turns))
    print("More importantly, this is a game of perfect information. For any given arrangement of coins, Player 1 can analyze the outcome of both of their possible first moves (taking the leftmost vs. taking the rightmost coin).")
    print("Player 1 can then choose the move that results in a better final score for them, forcing Player 2 to play in a less favorable scenario.")
    print("This strategic advantage means it is always better to be the first player.")

solve()
