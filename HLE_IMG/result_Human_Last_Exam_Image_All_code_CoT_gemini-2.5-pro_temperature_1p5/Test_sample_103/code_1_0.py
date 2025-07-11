def solve_chess_puzzle():
    """
    Analyzes the chess position to find the best move for White.
    """
    position_summary = {
        "white_pieces": "King h1, Rook a1, Queen d5, Knight f7, Pawns a2, b2, c3, g2, h2",
        "black_pieces": "King g8, Rooks d8 & f8, Queen h4, Bishop e4, Pawns a7, b7, g7, h7",
        "turn": "White",
        "main_threat": "Black's Queen on h4 threatens Qh1# (checkmate)."
    }

    moves = {
        "A": "Rf1", "B": "Qh4", "C": "Qxe4", "D": "h3",
        "E": "Ng5+", "F": "Qxd8", "G": "Ne5+", "H": "Nd6+",
        "I": "Nh6+", "J": "Nxd8+"
    }

    print("Analyzing the chess position to find the best move for White.")
    print("----------------------------------------------------------")
    print("Black's main threat is Qh1#, an immediate checkmate.")
    print("White must find a forcing move that wins the game while dealing with this threat.")
    print("\nLet's evaluate the options:\n")

    print("Moves like Qxd8, Nd6+, and Nxd8+ are blunders. They allow Black to play ...Qh1#, winning the game instantly.\n")
    print("Moves like Rf1 and h3 defend against the checkmate but are too passive. White has a winning attack and shouldn't defend.\n")
    print("Qxh4 trades queens, which is good for White, but there is a much stronger, more decisive continuation.\n")

    print("The best move initiates a forcing sequence that wins material.")
    print("The winning move is G. Ne5+.")
    print("\nHere is the winning combination:")
    print("1. Ne5+  (This check forces Black's King to move.)")
    print("   1... Kh8 (This is Black's only legal move.)")
    print("2. Ng6+  (A brilliant sacrifice to open the h-file for the White Queen!)")
    print("   2... hxg6 (Forced capture.)")
    print("3. Qxh4+ (White captures the Black Queen with check.)")
    print("   3... Kg8")
    print("\nAfter this sequence, White has won the black queen for a knight and has a completely winning position.")

solve_chess_puzzle()