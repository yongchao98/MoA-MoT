def find_winning_t_for_player2(limit):
    """
    Identifies and prints the total number of tokens T for which the second player
    has a winning strategy in the described token game.

    This occurs when T is a Fibonacci number (>= 2).
    """
    print(f"The second player has a winning strategy if the starting number of tokens T is a Fibonacci number.")
    print(f"The following are the winning values of T up to {limit}:")
    
    # We check T=2 as a base case, as the game is trivial.
    # Player 1 must take 1 (since 0 < x < 2), leaving 1. Player 2 takes the last token.
    print("T = 2 is a winning value for Player 2.")

    # We use the Fibonacci sequence defined by F(n) = F(n-1) + F(n-2).
    # We will start with the numbers that produce the sequence from 3 onwards.
    # f_prev corresponds to F(n-2), f_curr corresponds to F(n-1)
    f_prev, f_curr = 1, 2 

    while True:
        # Calculate the next Fibonacci number
        f_next = f_prev + f_curr
        
        if f_next > limit:
            break
            
        # The recurrence relation F_next = f_curr + f_prev serves as the 'equation'.
        print(f"T = {f_next} is a winning value for Player 2 (as in the Fibonacci sequence: {f_next} = {f_curr} + {f_prev})")
        
        # Update the values for the next iteration
        f_prev = f_curr
        f_curr = f_next

# We can find the winning values for T up to a certain limit, for example 500.
find_winning_t_for_player2(500)