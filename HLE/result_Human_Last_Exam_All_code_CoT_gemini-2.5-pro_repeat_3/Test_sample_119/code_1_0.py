import chess
import chess.pgn
import io

def solve_chess_mate():
    """
    This function analyzes a specific chess position to find the shortest
    mating sequence for black.
    """
    # The PGN of the game leading to the position
    pgn_string = "1. e4 c6 2. c4 e6 3. Nc3 Bb4 4. a3 Bxc3 5. bxc3 Nf6 6. f3 b6 7. d4 c5 8. d5 O-O 9. Ne2 Ba6 10. Nf4 exd5 11. cxd5 Bxf1 12. Rxf1 Re8 13. Kf2 d6 14. Kg1 Nfd7 15. c4 Ne5 16. Qc2 Qf6 17. g3 Nxc4 18. Rb1 Qd4+ 19. Kh1 Nd7 20. Rd1 Ne3 21. Bxe3 Qxe3 22. Rd3 Qxd3 23. Qxd3 Ne5 24. Qe2 a6 25. Rxb6 Reb8 26. Rxd6 Rb1+ 27. Kg2 Rab8 28. Kh3 R1b2 29. Qxa6 Nxf3 30. Qa5 Rxh2+ 31. Kg4 Ne5+ 32. Kf5 Re8 33. Rd8 g6+ 34. Kg5"

    # Read the PGN and set up the board to the final position
    pgn = io.StringIO(pgn_string)
    game = chess.pgn.read_game(pgn)
    board = game.end().board()

    # --- Find the shortest mate ---

    # Check for Mate in 1
    for move in board.legal_moves:
        board_copy = board.copy()
        board_copy.push(move)
        if board_copy.is_checkmate():
            # If a mate-in-1 is found, print it and we are done
            print(board.san(move))
            return

    # Check for Mate in 2
    for b1_move in board.legal_moves:
        board1 = board.copy()
        board1.push(b1_move)
        
        # We assume this move leads to a forced mate until proven otherwise
        is_forced_mate = True
        # Store one of the valid mating lines
        final_mate_line = []

        # White has no legal moves, means b1_move was mate. This is already
        # handled by the mate-in-1 check, but included for completeness.
        if not board1.legal_moves:
            continue

        # Iterate through all of White's possible replies
        for w1_move in board1.legal_moves:
            board2 = board1.copy()
            board2.push(w1_move)
            
            found_b2_mate = False
            # Check if Black has a mating response
            for b2_move in board2.legal_moves:
                board3 = board2.copy()
                board3.push(b2_move)
                if board3.is_checkmate():
                    # Found a mating move for this line.
                    final_mate_line = [board.san(b1_move), board1.san(w1_move), board2.san(b2_move)]
                    found_b2_mate = True
                    break # Stop searching for other mates in this sub-line
            
            if not found_b2_mate:
                # If for any of White's replies, Black cannot force mate,
                # then our initial move b1_move is not the solution.
                is_forced_mate = False
                break # Stop checking other White replies

        if is_forced_mate:
            # If we checked all of White's replies and found a mate for each,
            # we have found the shortest mating sequence.
            print(" ".join(final_mate_line))
            return
            
    print("No short mating sequence found.")

# Execute the solver
solve_chess_mate()