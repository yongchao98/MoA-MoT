def find_winning_t_for_p2():
    """
    This function determines and prints the initial number of tokens T for which the second player has a winning strategy.

    The analysis of this game reveals a connection to the Fibonacci numbers.
    A well-established result in combinatorial game theory (by A. Schwenk) states that
    the second player has a winning strategy if and only if the initial number of tokens, T,
    is a Fibonacci number.

    The relevant Fibonacci sequence starts with F_1=1, F_2=2.
    F_n = F_{n-1} + F_{n-2} for n > 2.
    Sequence: 1, 2, 3, 5, 8, 13, 21, ...
    """
    
    print("The second player has a winning strategy if and only if the total number of tokens T is a Fibonacci number.")
    print("The sequence of these winning numbers for the second player begins: 1, 2, 3, 5, 8, ...")
    
    limit = 1000
    print(f"\nThe values of T up to {limit} for which the second player wins are:")

    fib_numbers = []
    a, b = 1, 2
    while a <= limit:
        fib_numbers.append(a)
        a, b = b, a + b
    
    # Printing each number as requested
    for num in fib_numbers:
        print(num, end=' ')
    print()

find_winning_t_for_p2()
