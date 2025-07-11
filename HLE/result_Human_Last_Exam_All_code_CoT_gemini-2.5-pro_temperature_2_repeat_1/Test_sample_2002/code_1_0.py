def find_winning_t_for_player2(limit):
    """
    This function identifies the initial number of tokens T for which the
    second player has a guaranteed winning strategy.

    According to game theory on this specific variant of Nim, the second player
    has a winning strategy if and only if T is a Fibonacci number.
    The game requires T >= 2 for the first player to make a valid move.

    This code will print the Fibonacci numbers (F_n >= 2) up to a given limit,
    as these are the winning values for T. It also prints the recurrence
    relation for each number, as per the request to show the 'final equation'.
    """

    print("The second player has a winning strategy if and only if T is a Fibonacci number.")
    print(f"The winning values for T (up to {limit}) are:")

    # We use the standard Fibonacci sequence starting with F_1=1, F_2=1.
    # f_prev represents F_{n-2}, f_curr represents F_{n-1}
    f_prev = 1
    f_curr = 1

    # Loop to generate the next Fibonacci numbers
    while True:
        # F_n = F_{n-1} + F_{n-2}
        f_next = f_curr + f_prev

        if f_next > limit:
            break

        # T = f_next represents a winning starting position for Player 2.
        print(f"T = {f_next}, from the equation: {f_next} = {f_curr} + {f_prev}")

        # Update the values for the next iteration
        f_prev = f_curr
        f_curr = f_next

# Let's find the winning T values up to 200
find_winning_t_for_player2(200)
