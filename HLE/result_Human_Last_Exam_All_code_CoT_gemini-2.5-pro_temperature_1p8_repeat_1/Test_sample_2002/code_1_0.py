def find_winning_T_for_player2(limit):
    """
    This function finds the number of tokens T for which the second player
    has a winning strategy in the described token game.

    The winning strategy for the second player exists if and only if T
    is a Fibonacci number.

    Args:
        limit: The maximum value of T to check.
    """
    # The game is ill-defined for T < 1. If T=1, P1 has no move (0 < x < 1)
    # and loses, so P2 wins. Thus, T=1 is a winning number for P2.
    # The Fibonacci sequence that represents these numbers is 1, 2, 3, 5, ...
    
    winning_t_values = []
    # We initialize the Fibonacci sequence with a=1, b=2 to generate 1, 2, 3, 5,...
    # The value `a` will hold the current Fibonacci number in the loop.
    a, b = 1, 2
    
    while a <= limit:
        winning_t_values.append(a)
        # Move to the next Fibonacci number
        a, b = b, a + b
        
    print("The values of T for which the second player has a winning strategy are the Fibonacci numbers.")
    print(f"Here are the values up to {limit}:")
    print(", ".join(map(str, winning_t_values)))

# Set a limit for T
T_LIMIT = 1000
find_winning_T_for_player2(T_LIMIT)