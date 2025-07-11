def find_winning_positions_for_p2(limit):
    """
    This function identifies the initial number of tokens T for which the second player
    has a winning strategy in the described token game.

    The winning positions for the second player correspond to the Fibonacci numbers.
    This is a known result from the analysis of the game called Fibonacci Nim.

    Args:
    limit (int): The upper bound for the values of T to check.
    """

    # We determine the P2-winning numbers by experiment for small T: 1, 2, 3, 5, ...
    # This corresponds to the Fibonacci sequence starting with 1, 2.
    a, b = 1, 2
    
    # A list to store the values of T where the second player wins.
    p2_winning_t = []

    if a <= limit:
        p2_winning_t.append(a)

    while b <= limit:
        p2_winning_t.append(b)
        # Calculate the next Fibonacci number
        a, b = b, a + b
        
    print(f"The values of T for which the second player has a winning strategy are the Fibonacci numbers.")
    print(f"Here are these numbers up to {limit}:")
    
    # The prompt asks to output each number, so we format it as a clean, comma-separated string.
    output_string = ", ".join(map(str, p2_winning_t))
    print(output_string)

# Set a reasonable limit to display the winning numbers.
find_winning_positions_for_p2(200)