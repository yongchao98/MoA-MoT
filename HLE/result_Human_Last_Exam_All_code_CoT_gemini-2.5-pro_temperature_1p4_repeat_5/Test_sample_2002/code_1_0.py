def find_winning_t_for_p2(limit):
    """
    In the described token game (Fibonacci Nim), the second player has a
    winning strategy if and only if the initial number of tokens, T, is
    a Fibonacci number. This function prints these values up to a given limit.

    The sequence of losing positions starts with 1, 2, 3, 5, 8, ...
    """
    print("The values of T for which the second player has a winning strategy are the Fibonacci numbers.")
    print(f"Here are the winning values for T up to {limit}:")

    winning_values = []
    # Initialize with the first two Fibonacci numbers in the sequence for game theory
    a, b = 1, 2
    while a <= limit:
        winning_values.append(a)
        # Standard Fibonacci progression
        a, b = b, a + b

    # Print the numbers as requested
    output_str = ""
    for val in winning_values:
        output_str += str(val) + " "
    print(output_str.strip())


if __name__ == '__main__':
    # Set a reasonable limit for the number of tokens T
    T_LIMIT = 200
    find_winning_t_for_p2(T_LIMIT)