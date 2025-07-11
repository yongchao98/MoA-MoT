def solve_domino_game():
    """
    Analyzes the domino game to determine the winning player for different board sizes.

    The game is a misère play game equivalent to taking 1 or 2 items from a pile of size n.
    The last player to move loses.

    A position 'k' is a losing position for the current player if k mod 3 == 1.
    Otherwise, it is a winning position.

    The first player starts at position n. If n is a winning position, Player 1 wins.
    If n is a losing position, Player 2 wins.
    """

    print("Analyzing the game for different board sizes (n).")
    print("A position 'k' is a losing position for the current player if k mod 3 = 1.")
    print("----------------------------------------------------------------------")

    # List of board sizes to analyze
    board_sizes = [20, 24]
    results = {}

    for n in board_sizes:
        print(f"\nCase n = {n}:")
        
        # This is the final equation with numbers
        remainder = n % 3
        print(f"The final equation is: {n} % 3 = {remainder}")

        if remainder == 1:
            winner = "Second player"
            position_type = "Losing"
        else:
            winner = "First player"
            position_type = "Winning"
        
        results[n] = winner
        print(f"The remainder is {remainder}, so n={n} is a {position_type} position for the first player.")
        print(f"Therefore, the {winner} has a winning strategy.")

    print("\n----------------------------------------------------------------------")
    print("Summary:")
    print(f"When n = 20, the {results[20].lower()} has a winning strategy.")
    print(f"When n = 24, the {results[24].lower()} has a winning strategy.")

solve_domino_game()