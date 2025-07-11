def solve_spice_girls_chess_riddle():
    """
    This function solves the riddle by mapping the Spice Girls to a chessboard rank.
    """
    # Step 1: The order of the Spice Girls as referenced in the rap bridge.
    # "Em" -> Emma Bunton (Baby Spice)
    # "G" -> Geri Halliwell (Ginger Spice)
    # "V" -> Victoria Beckham (Posh Spice)
    # "me" -> Mel B, the rapper (Scary Spice)
    # "MC" -> Mel C, referenced in G's line (Sporty Spice)
    spice_girls_order = ["Baby", "Ginger", "Posh", "Scary", "Sporty"]

    # Step 2: The order of pieces on White's starting rank from queenside (a1) to kingside (h1).
    white_rank_pieces = ["Rook", "Knight", "Bishop", "Queen", "King", "Bishop", "Knight", "Rook"]

    # Step 3: Find the position of the White Queen. In a list, this is the index.
    # The pieces are on squares 1, 2, 3, 4, etc.
    # The Queen is on the 4th square.
    queen_position_index = white_rank_pieces.index("Queen")
    queen_position_number = queen_position_index + 1

    # Step 4: Find the Spice Girl at that position.
    target_member = spice_girls_order[queen_position_index]

    # Step 5: Print the logic and the final answer.
    # The prompt requires outputting the numbers in the "equation".
    print("The mapping of positions to Spice Girls is as follows:")
    for i, member in enumerate(spice_girls_order):
        position = i + 1
        piece = white_rank_pieces[i]
        print(f"Position {position} ({piece}'s square): {member} Spice")

    print(f"\nThe White Queen starts at position {queen_position_number}.")
    print(f"The member at position {queen_position_number} is {target_member} Spice.")
    print(f"\nTherefore, the answer is the word '{target_member}'.")


solve_spice_girls_chess_riddle()
<<<Scary>>>