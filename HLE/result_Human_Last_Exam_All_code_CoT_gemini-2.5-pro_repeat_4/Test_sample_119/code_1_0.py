import chess

def solve_chess_puzzle():
    """
    Solves the chess puzzle by finding the shortest mating sequence for Black.
    The code systematically checks for forced mates and prints the principal variation.
    """
    # FEN of the position after White's 34. Kg5. This is the starting point for our analysis.
    fen = "4r1k1/5p1p/6p1/Q1pPn1K1/4PN2/P5P1/7r/8 b - - 0 34"
    board = chess.Board(fen)

    # We will search for the shortest forced mate, starting with mate-in-1, then mate-in-2.
    # We prioritize checking moves as they are the most likely to lead to a quick mate.
    candidate_first_moves = []
    for move in board.legal_moves:
        board.push(move)
        is_check = board.is_check()
        board.pop()
        if is_check:
            # Insert checking moves at the beginning of the list to check them first.
            candidate_first_moves.insert(0, move)
        else:
            candidate_first_moves.append(move)

    # Iterate through Black's possible first moves to find one that forces mate.
    for black_move1 in candidate_first_moves:
        board.push(black_move1)
        
        # If this move is an immediate checkmate, we have found a mate in 1.
        if board.is_checkmate():
            board.pop()
            return [board.san(black_move1)]

        is_forced_mate = True
        main_line_white_move = None
        main_line_black_mate = None
        
        white_responses = list(board.legal_moves)
        if not white_responses:
            # This is stalemate, not checkmate.
            is_forced_mate = False

        # For the given Black move, check ALL of White's possible replies.
        for white_move in white_responses:
            board.push(white_move)
            
            found_mating_reply = False
            mating_reply = None
            # Check if Black has a mating move in this new position.
            for black_move2 in board.legal_moves:
                board.push(black_move2)
                if board.is_checkmate():
                    found_mating_reply = True
                    mating_reply = black_move2
                    board.pop() # backtrack black_move2
                    break
                board.pop() # backtrack black_move2
            
            board.pop() # backtrack white_move

            if not found_mating_reply:
                # If for any of White's moves, Black has no mating reply, then the initial
                # move does not force a mate.
                is_forced_mate = False
                break
            
            # To build the "main line" for the output, we favor White's move that
            # is a capture, as this is typically the most critical response to check.
            if board.is_capture(white_move) or main_line_white_move is None:
                main_line_white_move = white_move
                main_line_black_mate = mating_reply

        board.pop() # backtrack black_move1

        if is_forced_mate:
            # If we get here, it means black_move1 forces a mate in 2.
            # We will format this line for our answer.
            san_b1 = board.san(black_move1)
            san_w1 = board.san(main_line_white_move)
            
            # To get the correct notation for the final move (including check/mate symbols),
            # we replay the sequence on the board.
            board.push(black_move1)
            board.push(main_line_white_move)
            san_b2 = board.san(main_line_black_mate)
            
            return [san_b1, san_w1, san_b2]

    return None

# --- Main execution ---
mating_sequence = solve_chess_puzzle()

if mating_sequence:
    print(" ".join(mating_sequence))
else:
    # This part should not be reached for the given puzzle.
    print("No short mating sequence found.")
