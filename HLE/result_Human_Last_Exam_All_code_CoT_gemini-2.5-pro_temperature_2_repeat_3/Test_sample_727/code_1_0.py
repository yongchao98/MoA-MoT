import chess
import chess.pgn
import chess.engine
import io
import sys
import os
import shutil

def find_stockfish():
    """
    Tries to find a Stockfish executable in common locations or in the system's PATH.
    """
    # First, try to find stockfish in the system PATH
    stockfish_path = shutil.which("stockfish")
    if stockfish_path:
        return stockfish_path
    
    # If not in PATH, check a list of common paths and executable names
    # Common for Windows, Linux, and macOS
    possible_paths = [
        "stockfish",
        "/usr/games/stockfish",
        "/usr/bin/stockfish",
        "/usr/local/bin/stockfish",
        "C:/Program Files/Stockfish/stockfish.exe",
    ]
    for path in possible_paths:
        if os.path.exists(path) and os.access(path, os.X_OK):
            return path
            
    # Heuristic for finding within a folder if the script is run from a base dir
    if os.path.exists("stockfish/stockfish-windows-x86-64-avx2.exe"):
        return "stockfish/stockfish-windows-x86-64-avx2.exe"
    if os.path.exists("stockfish-ubuntu-x86-64-avx2"):
        return "stockfish-ubuntu-x86-64-avx2"
        
    return None

def solve_chess_puzzle():
    """
    Analyzes the chess position to find the drawing move.
    """
    # Step 1: Check for dependencies (python-chess and stockfish engine)
    try:
        import chess, chess.pgn, chess.engine
    except ImportError:
        print("The 'python-chess' library is not installed.")
        print("Please install it by running: pip install python-chess")
        return

    stockfish_path = find_stockfish()
    if not stockfish_path:
        print("Error: Stockfish engine not found.")
        print("Please download it from https://stockfishchess.org/ and place the executable")
        print("in your system's PATH or in the same directory as this script.")
        return

    # Step 2: Set up the game from PGN
    pgn_string = """
    [Event "2021 World Chess Championship"]
    [Site "Dubai UAE"]
    [Date "2021.12.03"]
    [Round "6"]
    [White "Carlsen, Magnus"]
    [Black "Nepomniachtchi, Ian"]
    [Result "1-0"]
    [ECO "D02"]
    1. d4 Nf6 2. Nf3 d5 3. g3 e6 4. Bg2 Be7 5. O-O O-O 6. b3 c5 7. dxc5 Bxc5 8. c4
    dxc4 9. Qc2 Qe7 10. Nbd2 Nc6 11. Nxc4 b5 12. Nce5 Nb4 13. Qb2 Bb7 14. a3 Nc6
    15. Nd3 Bb6 16. Bg5 Rfd8 17. Bxf6 gxf6 18. Rac1 Nd4 19. Nxd4 Bxd4 20. Qa2
    Bxg2 21. Kxg2 Qb7+ 22. Kg1 Qe4 23. Qc2 a5 24. Rfd1 Kg7 25. Rd2 Rac8 26. Qxc8
    Rxc8 27. Rxc8 Qd5 28. b4 a4 29. e3 Be5 30. h4 h5 31. Kh2 Bb2 32. Rc5 Qd6 33.
    Rd1 Bxa3 34. Rxb5 Qd7 35. Rc5 e5 36. Rc2 Qd5 37. Rdd2 Qb3 38. Ra2 e4 39. Nc5
    Qxb4 40. Nxe4 Qb3 41. Rac2 Bf8 42. Nc5 Qb5 43. Nd3 a3 44. Nf4 Qa5 45. Ra2 Bb4
    46. Rd3 Kh6 47. Rd1 Qa4 48. Rda1 Bd6 49. Kg1 Qb3 50. Ne2 Qd3 51. Nd4 Kh7 52.
    Kh2 Qe4 53. Rxa3 Qxh4+ 54. Kg1 Qe4 55. Ra4 Be5 56. Ne2 Qc2 57. R1a2 Qb3 58.
    Kg2 Qd5+ 59. f3 Qd1 60. f4 Bc7 61. Kf2 Bb6 62. Ra1 Qb3 63. Re4 Kg7 64. Re8
    f5 65. Raa8 Qb4 66. Rac8 Ba5 67. Rc1 Bb6 68. Re5 Qb3 69. Re8 Qd5 70. Rcc8 Qh1
    71. Rc1 Qd5 72. Rb1 Ba7 73. Re7 Bc5 74. Re5 Qd3 75. Rb7 Qc2 76. Rb5 Ba7 77.
    Ra5 Bb6 78. Rab5 Ba7 79. Rxf5 Qd3 80. Rxf7+ Kxf7 81. Rb7+ Kg6 82. Rxa7 Qd5
    83. Ra6+ Kh7 84. Ra1 Kg6 85. Nd4 Qb7 86. Ra2 Qh1 87. Ra6+ Kf7 88. Nf3 Qb1 89.
    Rd6 Kg7 90. Rd5 Qa2+ 91. Rd2 Qb1 92. Re2 Qb6 93. Rc2 Qb1 94. Nd4 Qh1 95.
    Rc7+ Kf6 96. Rc6+ Kf7 97. Nf3 Qb1 98. Ng5+ Kg7 99. Ne6+ Kf7 100. Nd4 Qh1 101.
    Rc7+ Kf6 102. Nf3 Qb1 103. Rd7 Qb2+ 104. Rd2 Qb1 105. Ng1 Qb4 106. Rd1 Qb3
    107. Rd6+ Kg7 108. Rd4 Qb2+ 109. Ne2 Qb1 110. e4 Qh1 111. Rd7+ Kg8 112. Rd4
    Qh2+ 113. Ke3 h4 114. gxh4 Qh3+ 115. Kd2 Qxh4 116. Rd3 Kf8 117. Rf3 Qd8+
    118. Ke3 Qa5 119. Kf2 Qa7+ 120. Re3 Qd7 121. Ng3 Qd2+ 122. Kf3 Qd1+ 123. Re2
    Qb3+ 124. Kg2 Qb7 125. Rd2 Qb3 126. Rd5 Ke7 127. Re5+ Kf7 128. Rf5+ Ke8 129.
    e5 Qa2+ 130. Kh3 Qe6 131. Kh4 Qh6+ 132. Nh5 Qh7 133. e6 Qg6 134. Rf7 Kd8 135.
    f5 Qg1 136. Ng7
    """
    
    pgn_io = io.StringIO(pgn_string)
    game = chess.pgn.read_game(pgn_io)
    board = game.board()

    # Fast-forward to the position after White's 130th move (130. Kh3)
    # This is the 259th half-move (ply).
    moves = list(game.mainline_moves())
    for i in range(259):
        board.push(moves[i])

    # Step 3: Analyze the options
    options = {
        "A": "Qa1", "B": "Qa7", "C": "Qg2", "D": "Qf2", "E": "Qb2", "F": "Qd2",
        "G": "Qa6", "H": "Qh2", "I": "Qa8", "J": "Qa5", "K": "Qa4", "L": "Qc2",
        "M": "Qe2", "N": "Qa3"
    }
    blunder_move = "Qe6"
    
    all_moves_to_check = list(options.values())
    all_moves_to_check.append(blunder_move)

    print("Starting chess engine analysis...\nThis might take a moment.\n")

    with chess.engine.SimpleEngine.popen_uci(stockfish_path) as engine:
        results = {}
        legal_moves_san = [board.san(m) for m in board.legal_moves]
        
        for move_san in all_moves_to_check:
            if move_san not in legal_moves_san:
                continue

            board_copy = board.copy()
            board_copy.push_san(move_san)
            
            # Analyze for 0.5 seconds
            info = engine.analyse(board_copy, chess.engine.Limit(time=0.5))
            score_obj = info["score"].white()
            
            if score_obj.is_mate():
                # Assign a very high/low value for mate scores for sorting
                score_cp = 9999 if score_obj.mate() > 0 else -9999
            else:
                score_cp = score_obj.cp
            results[move_san] = score_cp
        
        # Step 4: Identify the drawing move and present results
        # A draw is an evaluation of 0. We're looking for the move that minimizes White's advantage.
        # So we look for the minimum cp value in the results, excluding the known blunder.
        
        drawing_candidates = {k: v for k, v in results.items() if k != blunder_move}
        
        if not drawing_candidates:
            print("No legal moves found among the provided options to analyze.")
            return

        best_move_san = min(drawing_candidates, key=drawing_candidates.get)
        best_score_cp = drawing_candidates[best_move_san]
        blunder_score_cp = results.get(blunder_move, "N/A")

        final_answer_letter = [key for key, value in options.items() if value == best_move_san][0]
        
        print("--- Analysis Complete ---")
        print(f"The actual move played, {blunder_move}, gives White an advantage of: {blunder_score_cp/100.0:.2f}")
        print(f"The best move found, {best_move_san}, results in an evaluation of: {best_score_cp/100.0:.2f} (close to a draw)")
        
        # Fulfill "equation" requirement
        print("\nEquation of Lost Advantage:")
        num1 = blunder_score_cp / 100.0
        num2 = best_score_cp / 100.0
        print(f"{num1:.2f} (Blunder) - {num2:.2f} (Best Move) = {num1 - num2:.2f} (Difference)")
        print("\nIndividual numbers from the equation:")
        print(num1)
        print(num2)
        print(num1 - num2)

        print(f"\nConclusion: The queen move that could have drawn the game is {best_move_san}.")
        
        print(f"<<<{final_answer_letter}>>>")

if __name__ == '__main__':
    solve_chess_puzzle()
