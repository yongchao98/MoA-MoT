def find_winner(n):
    """
    Determines the winner of the domino game for a board of size n.
    The last player to move loses.

    Args:
      n: The size of the n x 1 board.

    Returns:
      A string indicating which player has a winning strategy.
    """
    print(f"Analyzing the game for n = {n}:")

    # The winning/losing state is determined by the remainder of n divided by 3.
    # The 'final equation' is this modulus operation.
    remainder = n % 3
    print(f"{n} % 3 = {remainder}")

    # A position 'n' is a losing position for the current player if and only if n % 3 == 1.
    # In this case, the second player has a winning strategy.
    # Otherwise, it's a winning position for the first player.
    if remainder == 1:
        winner = "the second player"
    else:
        winner = "the first player"

    print(f"Result: When n = {n}, {winner} has a winning strategy.\n")
    return winner

if __name__ == '__main__':
    find_winner(20)
    find_winner(24)