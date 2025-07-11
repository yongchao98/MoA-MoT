def find_winning_t_for_player2(limit=100):
    """
    This function identifies and prints the initial number of tokens, T,
    for which the second player has a guaranteed winning strategy in the game of Fibonacci Nim.

    The winning strategy for the second player exists if and only if T is a Fibonacci number.
    The game requires the first player to take a positive number of tokens less than T,
    so the game cannot start with T=1. Thus, we look for Fibonacci numbers greater than 1.
    """
    print(f"The second player has a winning strategy when T is a Fibonacci number greater than 1.")
    print(f"The winning values for T (up to {limit}) are shown below as sums of the previous two Fibonacci numbers:")

    # Initialize with the first two Fibonacci numbers to generate the sequence.
    # F_1 = 1, F_2 = 1
    f_prev = 1  # Represents F_{n-2}
    f_curr = 1  # Represents F_{n-1}

    while True:
        # Calculate the next Fibonacci number
        t_win = f_curr + f_prev
        
        if t_win > limit:
            break

        # The prompt requires showing the equation for each number
        print(f"T = {t_win} = {f_curr} + {f_prev}")

        # Update the values for the next iteration
        f_prev = f_curr
        f_curr = t_win

if __name__ == '__main__':
    # Set an upper limit for T to demonstrate the solution
    find_winning_t_for_player2(limit=200)
