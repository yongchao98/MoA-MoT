def find_winning_t_for_player2():
    """
    This function determines and prints the initial number of tokens T
    for which the second player has a winning strategy in the described game.

    The game is a variant of Nim known as Fibonacci Nim. The winning strategy
    is determined by the Fibonacci sequence. The second player has a winning
    strategy if and only if T is a Fibonacci number. The relevant sequence
    starts with 1 and 2.

    This code generates and prints these Fibonacci numbers up to a limit.
    """

    print("The second player has a winning strategy if and only if T is a Fibonacci number.")
    print("The relevant Fibonacci sequence starts with F_1=1, F_2=2.")
    print("The winning values for T are:")

    limit = 1000  # We'll list the winning numbers up to this limit.
    a, b = 1, 2

    if a <= limit:
        print(a)
    if b <= limit:
        print(b)

    c = a + b
    while c <= limit:
        # The prompt asks to "output each number in the final equation"
        # which we interpret as showing the recurrence relation.
        print(f"{c} (since {c} = {b} + {a})")
        a, b = b, c
        c = a + b

find_winning_t_for_player2()