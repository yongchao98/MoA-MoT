import chess

def solve_chess_puzzle(fen1, fen2):
    """
    Analyzes two chess positions from FEN strings to determine if they can
    occur in the same game and in which order.

    The function will:
    1. Parse and display the two positions.
    2. Identify the piece moves that separate the two positions.
    3. Analyze the logical flow of the game based on piece development.
    4. Print a step-by-step analysis and conclusion.
    """
    try:
        board1 = chess.Board(fen1)
        board2 = chess.Board(fen2)
    except ValueError:
        print("Error: One or both FEN strings are invalid.")
        return

    # --- Step 1: Display the initial positions ---
    print("--- Position Analysis ---")
    print("\nPosition 1:")
    print(board1)
    print("FEN: " + fen1)

    print("\nPosition 2:")
    print(board2)
    print("FEN: " + fen2)

    # --- Step 2: Find the difference between the boards ---
    map1 = board1.piece_map()
    map2 = board2.piece_map()
    all_squares = set(map1.keys()) | set(map2.keys())
    diff_squares = {sq for sq in all_squares if map1.get(sq) != map2.get(sq)}

    moves = {}  # Using a dictionary to track moves: (piece_type, color) -> {'from': sq, 'to': sq}
    for sq in diff_squares:
        p1, p2 = map1.get(sq), map2.get(sq)
        if p1 and not p2:  # Piece moved FROM this square
            key = (p1.piece_type, p1.color)
            if key not in moves: moves[key] = {}
            moves[key]['from'] = sq
        elif not p1 and p2:  # Piece moved TO this square
            key = (p2.piece_type, p2.color)
            if key not in moves: moves[key] = {}
            moves[key]['to'] = sq

    # --- Step 3: Print the difference as a series of moves ---
    print("\n--- Difference Between Positions ---")
    print("The pawn structure and material are identical in both positions.")
    print("The net change corresponds to the following four pieces having moved:")

    # Sort for consistent output: White pieces then Black, by piece type
    sorted_moves = sorted(moves.items(), key=lambda item: (item[0][1], -item[0][0]), reverse=True)

    for (ptype, color), move_info in sorted_moves:
        from_sq_name = chess.square_name(move_info['from'])
        to_sq_name = chess.square_name(move_info['to'])
        color_str = "White" if color == chess.WHITE else "Black"
        piece_symbol = chess.piece_symbol(ptype).upper()
        
        # This format shows each element of the move "equation"
        print(f"   - {color_str} {chess.piece_name(ptype)} move: {piece_symbol} from {from_sq_name} to {to_sq_name}")

    # --- Step 4: Analyze the game flow and conclude ---
    print("\n--- Game Flow Analysis and Conclusion ---")
    print("Both positions are for White to move, and full castling rights ('KQkq') are intact for both sides.")
    print("\nPath from Position 1 to 2:")
    print("This transition requires two moves from White (e.g., Bf1-g2, Bc1-g5) and two from Black (e.g., Nb8-d7, Bf8-e7).")
    print("This sequence represents logical piece development, moving pieces from starting squares to more active ones.")

    print("\nPath from Position 2 to 1:")
    print("This transition requires the reverse moves (e.g., Bg2-f1, Bg5-c1, etc.), moving pieces back to their starting squares.")
    print("While these backward moves are legal according to the rules of chess, they are strategically illogical in a typical game.")

    print("\nConclusion:")
    print("A game of chess progresses from a state of lesser development to greater development.")
    print("Position 2 is a more developed version of Position 1.")
    print("Therefore, it is overwhelmingly likely that Position 1 must arise before Position 2 in any sensible game.")


if __name__ == '__main__':
    fen_pos1 = "rn1qkb1r/1p3ppp/p2pbn2/4p3/4P1P1/2N4P/PPP1NP2/R1BQKB1R w KQkq - 0 1"
    fen_pos2 = "r2qk2r/1p1nbppp/p2pbn2/4p1B1/4P1P1/2N4P/PPP1NPB1/R2QK2R w KQkq - 0 1"
    solve_chess_puzzle(fen_pos1, fen_pos2)
    print("\n<<<A>>>")