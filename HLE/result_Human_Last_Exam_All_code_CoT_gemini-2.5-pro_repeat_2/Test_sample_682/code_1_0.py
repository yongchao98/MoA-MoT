def solve_domino_game():
    """
    Determines the winner of a misere domino game on an n x 1 board.

    The game is impartial, and the last player to move loses. The winning
    strategy depends on the size of the board, n. A position is losing
    for the current player if n % 3 == 1. Otherwise, it is a winning position.
    """

    print("Analyzing the domino game with misere play (last player loses)...")
    print("A starting board of size 'n' is a losing position for the first player if and only if n mod 3 = 1.")
    print("-" * 40)

    # Case n = 20
    n1 = 20
    print(f"Analyzing for n = {n1}:")
    remainder1 = n1 % 3
    # This prints the equation as requested
    print(f"The calculation is: {n1} % 3 = {remainder1}")

    if remainder1 == 1:
        winner1 = "the second player"
    else:
        winner1 = "the first player"
    print(f"Result: When n = {n1}, {winner1} has a winning strategy.")
    print("-" * 40)

    # Case n = 24
    n2 = 24
    print(f"Analyzing for n = {n2}:")
    remainder2 = n2 % 3
    # This prints the equation as requested
    print(f"The calculation is: {n2} % 3 = {remainder2}")

    if remainder2 == 1:
        winner2 = "the second player"
    else:
        winner2 = "the first player"
    print(f"Result: When n = {n2}, {winner2} has a winning strategy.")
    print("-" * 40)

    # Final Conclusion
    print("Summary of results:")
    print(f"For n = {n1}: {winner1} wins.")
    print(f"For n = {n2}: {winner2} wins.")
    print("\nThis corresponds to the statement: 'When n = 20, the first player has a winning strategy. When n = 24, the first player has a winning strategy.'")

# Execute the analysis
solve_domino_game()