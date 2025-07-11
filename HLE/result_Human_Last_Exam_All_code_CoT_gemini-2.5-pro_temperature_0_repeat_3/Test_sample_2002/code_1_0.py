def find_winning_t_for_p2(limit):
    """
    This function finds the values of T for which the second player has a winning strategy.
    According to game theory, these values are the Fibonacci numbers.

    The problem asks to show the equation for each number, so we will print
    F_n = F_{n-1} + F_{n-2} for each Fibonacci number in the sequence.
    The sequence of losing positions for Player 1 is 1, 2, 3, 5, 8, ...
    This corresponds to the Fibonacci sequence F_n for n >= 2 (if F_0=0, F_1=1).
    """
    print(f"The values of T up to {limit} for which the second player has a winning strategy are the Fibonacci numbers:")
    print("-" * 80)

    # We start with F_1=1, F_2=1 to generate the sequence.
    # The first losing T is F_2=1.
    # The second losing T is F_3=2.
    # The third losing T is F_4=3.
    # And so on.

    # Let's use f_prev for F_{n-2} and f_curr for F_{n-1}
    f_prev, f_curr = 1, 1  # Start with F_1, F_2

    # Print the first losing number, F_2=1
    if f_curr <= limit:
        print(f_curr)
    else:
        return

    # Loop to generate and print the rest of the Fibonacci numbers
    while True:
        f_next = f_curr + f_prev
        if f_next > limit:
            break
        print(f"{f_next} = {f_curr} + {f_prev}")
        f_prev, f_curr = f_curr, f_next

# Set a limit for T and run the function
find_winning_t_for_p2(1000)