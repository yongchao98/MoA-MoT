def solve_domino_game(n):
    """
    Determines the winner of the domino game for a board of size n.
    The last player to play loses.
    A position n is losing for the current player if n % 3 == 1.
    Otherwise, it is a winning position.
    """
    print(f"For n = {n}:")
    
    remainder = n % 3
    
    # Print the equation for the modulo operation
    # The format f"{n} % 3 = {remainder}" shows the equation and its result.
    print(f"{n} % 3 = {remainder}")
    
    if remainder == 1:
        # This is a losing position for the first player.
        print("The second player has a winning strategy.\n")
    else:
        # This is a winning position for the first player.
        print("The first player has a winning strategy.\n")

# Analyze the cases from the problem
solve_domino_game(20)
solve_domino_game(24)
