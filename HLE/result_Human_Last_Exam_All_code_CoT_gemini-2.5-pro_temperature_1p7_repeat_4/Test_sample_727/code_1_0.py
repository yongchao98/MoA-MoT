def find_drawing_move():
    """
    This function analyzes a critical position from game 6 of the 2021 World Chess Championship
    and identifies the correct drawing move for Black.
    """
    print("Analyzing the position before Black's move 130.")
    print("--------------------------------------------------")
    print("The position is as follows:")
    print("  - White: King on h3, Rook on f5, Knight on g3, Pawn on e5.")
    print("  - Black: King on e8, Queen on a2.")
    print("It is Black's turn to move.\n")

    print("In the game, Black played 130... Qe6. This was the losing move because it allows the White King to escape.")
    print("For example: 130... Qe6? 131. Kh4! The White King is now safer, and the Knight can join the defense via Nh5. White's pawn will eventually promote.\n")

    print("The correct move for Black is one that forces a perpetual check. Let's examine the drawing move, C. Qg2.\n")
    
    print("The move 130... Qg2+ initiates a forced series of checks from which the White King cannot escape. This can be seen as the 'final equation' for a draw:\n")
    
    # Displaying the key drawing line with move numbers
    print("The drawing sequence:")
    print("Move 130: Black plays ...Qg2+")
    print("Move 131: White must play Kh4. (Note: 131. Nh5?? would be checkmate after ...Qxh5#)")
    print("Move 131: Black responds with ...Qh2+")
    print("Move 132: White is forced to play Kg4.")
    print("Move 132: Black continues with ...Qg2+, forcing the White king back to h4.")
    
    print("\nThis pattern of ...Qh2+ and ...Qg2+ forces a draw by threefold repetition, as White has no way to avoid the checks and make progress.")
    print("Therefore, the move that could have drawn the game is Qg2.")

find_drawing_move()
print("\n<<<C>>>")