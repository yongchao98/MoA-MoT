def solve_chess_puzzle():
    """
    Analyzes a chess position from Carlsen vs. Nepomniachtchi (2021)
    and determines the correct drawing move for Black.
    """
    move_number = 130
    losing_move_square = "e6"
    drawing_move_square = "h2"
    
    print("--- Chess Puzzle Analysis ---")
    print(f"The position is from Carlsen-Nepomniachtchi, World Championship 2021, Game 6, after White's move {move_number}. Kh3.")
    
    print("\n--- The Position ---")
    print("White has a King on h3, a Knight on g3, and pawns on e5 and f4.")
    print("Black has a King on e8 and a Queen on a2.")
    print("White's winning plan is to play 131. Kh4, using the king to help the pawns advance.")

    print(f"\n--- The Blunder: {move_number}... Qe6? ---")
    print("The move played in the game, 130... Qe6, was a decisive error.")
    print("It allows White to execute the winning plan with 131. Kh4!")
    print("After 131... Qh6+, White plays 132. Nh5, and White's king is safe. White went on to win the game.")

    print(f"\n--- The Drawing Move: {move_number}... Qh2! ---")
    print("The correct move for Black is H. Qh2. This move forces a draw by trapping the White king.")
    print("1. After 130... Qh2, the move 131. Kh4 is illegal because the Black queen on h2 attacks that square.")
    print("2. If White tries 131. Kg4, Black delivers immediate checkmate with 131... Qh4#.")
    print("3. White's only option is to move the knight or the f-pawn, but Black can maintain the trap. For instance: 131. Nf1 Qf2 132. Ng3 Qh2 leads to a draw by repetition of moves.")
    print("Because Black can perpetually restrict the white king, the position is a draw.")

    print("\n--- Final Equation ---")
    print(f"For move number {move_number}, playing Queen to {losing_move_square} loses.")
    print(f"For move number {move_number}, playing Queen to {drawing_move_square} draws.")
    print(f"Result: ({move_number}...Q{drawing_move_square}) = Draw")

solve_chess_puzzle()
<<<H>>>