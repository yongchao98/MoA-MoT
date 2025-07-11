def find_winning_t_for_p2(limit):
    """
    This function finds the values of T for which the second player has a winning strategy
    in the described token game. This occurs when T is a Fibonacci number.
    
    Args:
    limit (int): The upper bound for the number of tokens T.
    """
    print(f"The values of T up to {limit} for which the second player has a winning strategy are:")
    # The winning T values are the Fibonacci numbers (1, 2, 3, 5, 8, ...).
    # We start with a=1, b=2 to generate the sequence used in the game theory proof (1,2,3,5,...),
    # but we will print T=1 explicitly as it's a winning position for P2 as well.
    print(1)
    a, b = 2, 3
    while a <= limit:
        print(a)
        a, b = b, a + b

if __name__ == '__main__':
    # We will find the winning T values up to 100.
    limit = 100
    find_winning_t_for_p2(limit)