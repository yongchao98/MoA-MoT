import chess
import chess.pgn
import io

def solve_chess_puzzle():
    """
    This function sets up a chess board from a given PGN,
    finds the shortest mating sequence for the current player,
    and prints the result.
    """
    # The PGN string of the game moves
    pgn_string = "1. e4 c6 2. c4 e6 3. Nc3 Bb4 4. a3 Bxc3 5. bxc3 Nf6 6. f3 b6 7. d4 c5 8. d5 O-O 9. Ne2 Ba6 10. Nf4 exd5 11. cxd5 Bxf1 12. Rxf1 Re8 13. Kf2 d6 14. Kg1 Nfd7 15. c4 Ne5 16. Qc2 Qf6 17. g3 Nxc4 18. Rb1 Qd4+ 19. Kh1 Nd7 20. Rd1 Ne3 21. Bxe3 Qxe3 22. Rd3 Qxd3 23. Qxd3 Ne5 24. Qe2 a6 25. Rxb6 Reb8 26. Rxd6 Rb1+ 27. Kg2 Rab8 28. Kh3 R1b2 29. Qxa6 Nxf3 30. Qa5 Rxh2+ 31. Kg4 Ne5+ 32. Kf5 Re8 33. Rd8 g6+ 34. Kg5"

    # Use python-chess to read the PGN and set up the board
    pgn = io.StringIO(pgn_string)
    game = chess.pgn.read_game(pgn)
    board = game.end().board()

    # The problem asks for the shortest mating sequence. We'll look for a mate in 1.
    mating_sequence = []

    # Iterate through all legal moves for Black
    for move in board.legal_moves:
        # Create a temporary copy of the board to test the move
        temp_board = board.copy()
        temp_board.push(move)

        # Check if the move results in a checkmate
        if temp_board.is_checkmate():
            # If it is a checkmate, we've found our move.
            # Get the move in Standard Algebraic Notation (e.g., "Kg7#")
            mating_move_san = board.san(move)
            mating_sequence.append(mating_move_san)
            break  # Exit the loop as we've found the shortest mate

    # Print the result
    if mating_sequence:
        print("The shortest mating sequence is:")
        # The result is a single move for a mate-in-1
        print(' '.join(mating_sequence))
    else:
        print("No mate-in-one found.")

solve_chess_puzzle()