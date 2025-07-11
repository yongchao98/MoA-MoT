import chess
import chess.pgn
import io

def find_shortest_mate():
    """
    Analyzes the given chess game to find the shortest mating sequence for black.
    """
    # The PGN moves of the game provided by the user.
    pgn_string = """
    1. e4 c6 2. c4 e6 3. Nc3 Bb4 4. a3 Bxc3 5. bxc3 Nf6 6. f3 b6 7. d4 c5
    8. d5 O-O 9. Ne2 Ba6 10. Nf4 exd5 11. cxd5 Bxf1 12. Rxf1 Re8 13. Kf2 d6
    14. Kg1 Nfd7 15. c4 Ne5 16. Qc2 Qf6 17. g3 Nxc4 18. Rb1 Qd4+ 19. Kh1 Nd7
    20. Rd1 Ne3 21. Bxe3 Qxe3 22. Rd3 Qxd3 23. Qxd3 Ne5 24. Qe2 a6 25. Rxb6 Reb8
    26. Rxd6 Rb1+ 27. Kg2 Rab8 28. Kh3 R1b2 29. Qxa6 Nxf3 30. Qa5 Rxh2+
    31. Kg4 Ne5+ 32. Kf5 Re8 33. Rd8 g6+ 34. Kg5
    """

    # Use the python-chess library to load the game and get the final board state.
    pgn = io.StringIO(pgn_string)
    game = chess.pgn.read_game(pgn)
    board = game.end().board()

    # --- Search for Mate in 1 ---
    for move in board.legal_moves:
        if board.is_legal(move):
            board.push(move)
            if board.is_checkmate():
                # Found a mate in 1.
                print(board.san(move))
                return
            board.pop()

    # --- Search for Mate in 2 ---
    solution_path = None
    for black_move1 in board.legal_moves:
        # Get the notation for the first move
        san1 = board.san(black_move1)
        board.push(black_move1)
        
        is_forced_mate = True
        white_responses = list(board.legal_moves)
        
        # If White has no legal moves, it's already mate (or stalemate)
        if not white_responses:
            is_forced_mate = False # Should have been caught by mate-in-1 check
        
        path_for_this_branch = []

        for white_move in white_responses:
            san2 = board.san(white_move)
            board.push(white_move)
            
            found_mating_reply = False
            for black_move2 in board.legal_moves:
                # Check if this move is checkmate
                board.push(black_move2)
                if board.is_checkmate():
                    san3 = board.san(black_move2)
                    found_mating_reply = True
                    # This is a valid mating line. Since we only need one, we can store it.
                    # For a forced mate, all white replies must have a mating response.
                    path_for_this_branch = [san1, san2, san3]
                    board.pop()
                    break
                board.pop()
            
            board.pop() # Undo white_move
            
            if not found_mating_reply:
                # This white move escapes the mate, so black_move1 is not a forced mate.
                is_forced_mate = False
                break
        
        board.pop() # Undo black_move1
        
        if is_forced_mate and path_for_this_branch:
            solution_path = path_for_this_branch
            break # Found the shortest forced mate sequence

    if solution_path:
        for move_san in solution_path:
            print(move_san)
    else:
        print("No short mating sequence found.")

# Execute the function to find and print the mate.
find_shortest_mate()
<<<f6+
Kh6
R8e6#>>>