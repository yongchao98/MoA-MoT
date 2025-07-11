import chess

def solve_chess_puzzle():
    """
    This function sets up the described chess position and verifies the
    shortest mating sequence for White.
    """
    # Step 1: Define the position using the FEN string derived from the user's description.
    # The FEN represents the placement of all pieces, whose turn it is (White 'w'),
    # castling rights (White can castle King and Queen side 'KQ', Black cannot '-'),
    # and move counters.
    fen = "rn3rk1/pbpqp1pp/1p1p1b2/4N2Q/3N1P2/3B4/PP1P1PP1/R3K2R w KQ - 0 1"
    board = chess.Board(fen)

    print("--- Chess Puzzle Analysis ---")
    print("Initial Board Position:")
    print(board)
    print("\nWhite to move.")

    # Step 2: Define the moves for the discovered mating sequence.
    # The moves are given in Standard Algebraic Notation (SAN).
    # This is a mate in 2 moves.
    number_of_moves = 2
    move1_white_san = "Ng6+"
    move1_black_san = "Kh8"
    move2_white_san = "Qxh7#"

    # Step 3: Programmatically verify the sequence.
    # We use a copy of the board to play through the moves.
    board_copy = board.copy()

    # Apply White's first move
    move1_w = board_copy.parse_san(move1_white_san)
    board_copy.push(move1_w)

    # Check that Black's reply is forced (only one legal move)
    is_forced = False
    if len(list(board_copy.legal_moves)) == 1:
        # The only legal move for black is Kh8
        is_forced = True
        board_copy.push(list(board_copy.legal_moves)[0])

    # Apply White's second move (the checkmate)
    move2_w = board_copy.parse_san(move2_white_san)
    board_copy.push(move2_w)
    
    # Step 4: Confirm that the final position is checkmate and print the results.
    if board_copy.is_checkmate() and is_forced:
        print(f"\nResult: White can deliver a forced checkmate in {number_of_moves} moves.")
        print("\nThe final equation for the checkmate is:")
        # The prompt asks to "output each number in the final equation".
        # The following lines print the numbered move sequence.
        print(f"1. {move1_white_san} {move1_black_san}")
        print(f"2. {move2_white_san}")
    else:
        print("\nCould not confirm the checkmate sequence. Please re-evaluate the position.")

if __name__ == '__main__':
    # To run this code, you need to install the python-chess library:
    # pip install python-chess
    solve_chess_puzzle()
