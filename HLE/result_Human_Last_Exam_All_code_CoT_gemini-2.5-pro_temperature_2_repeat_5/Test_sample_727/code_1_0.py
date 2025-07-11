import chess
import chess.pgn
import chess.engine
import io
import re

def find_drawing_move():
    """
    This function analyzes a specific position from a chess game to find a drawing move.
    It requires the 'python-chess' library and a working Stockfish engine installation.
    """
    # The PGN of the game as provided in the problem.
    pgn_string = """
    1.d4 Nf6 2. Nf3 d5 3. g3 e6 4. Bg2 Be7 5. O-O O-O 6. b3 c5 7. dxc5 Bxc5 8. c4 dxc4 9. Qc2 Qe7 10. Nbd2 Nc6 11. Nxc4 b5 12. Nce5 Nb4 13. Qb2 Bb7 14. a3 Nc6 15. Nd3 Bb6 16. Bg5 Rfd8 17. Bxf6 gxf6 18. Rac1 Nd4 19. Nxd4 Bxd4 20. Qa2 Bxg2 21. Kxg2 Qb7+ 22. Kg1 Qe4 23. Qc2 a5 24. Rfd1 Kg7 25. Rd2 Rac8 26. Qxc8 Rxc8 27. Rxc8 Qd5 28. b4 a4 29. e3 Be5 30. h4 h5 31. Kh2 Bb2 32. Rc5 Qd6 33. Rd1 Bxa3 34. Rxb5 Qd7 35. Rc5 e5 36. Rc2 Qd5 37. Rdd2 Qb3 38. Ra2 e4 39. Nc5 Qxb4 40. Nxe4 Qb3 41. Rac2 Bf8 42. Nc5 Qb5 43. Nd3 a3 44. Nf4 Qa5 45. Ra2 Bb4 46. Rd3 Kh6 47. Rd1 Qa4 48. Rda1 Bd6 49. Kg1 Qb3 50. Ne2 Qd3 51. Nd4 Kh7 52. Kh2 Qe4 53. Rxa3 Qxh4+ 54. Kg1 Qe4 55. Ra4 Be5 56. Ne2 Qc2 57. R1a2 Qb3 58. Kg2 Qd5+ 59. f3 Qd1 60. f4 Bc7 61. Kf2 Bb6 62. Ra1 Qb3 63. Re4 Kg7 64. Re8 f5 65. Raa8 Qb4 66. Rac8 Ba5 67. Rc1 Bb6 68. Re5 Qb3 69. Re8 Qd5 70. Rcc8 Qh1 71. Rc1 Qd5 72. Rb1 Ba7 73. Re7 Bc5 74. Re5 Qd3 75. Rb7 Qc2 76. Rb5 Ba7 77. Ra5 Bb6 78. Rab5 Ba7 79. Rxf5 Qd3 80. Rxf7+ Kxf7 81. Rb7+ Kg6 82. Rxa7 Qd5 83. Ra6+ Kh7 84. Ra1 Kg6 85. Nd4 Qb7 86. Ra2 Qh1 87. Ra6+ Kf7 88. Nf3 Qb1 89. Rd6 Kg7 90. Rd5 Qa2+ 91. Rd2 Qb1 92. Re2 Qb6 93. Rc2 Qb1 94. Nd4 Qh1 95. Rc7+ Kf6 96. Rc6+ Kf7 97. Nf3 Qb1 98. Ng5+ Kg7 99. Ne6+ Kf7 100. Nd4 Qh1 101. Rc7+ Kf6 102. Nf3 Qb1 103. Rd7 Qb2+ 104. Rd2 Qb1 105. Ng1 Qb4 106. Rd1 Qb3 107. Rd6+ Kg7 108. Rd4 Qb2+ 109. Ne2 Qb1 110. e4 Qh1 111. Rd7+ Kg8 112. Rd4 Qh2+ 113. Ke3 h4 114. gxh4 Qh3+ 115. Kd2 Qxh4 116. Rd3 Kf8 117. Rf3 Qd8+ 118. Ke3 Qa5 119. Kf2 Qa7+ 120. Re3 Qd7 121. Ng3 Qd2+ 122. Kf3 Qd1+ 123. Re2 Qb3+ 124. Kg2 Qb7 125. Rd2 Qb3 126. Rd5 Ke7 127. Re5+ Kf7 128. Rf5+ Ke8 129. e5 Qa2+ 130. Kh3 Qe6 131. Kh4 Qh6+ 132. Nh5 Qh7 133. e6 Qg6 134. Rf7 Kd8 135. f5 Qg1 136. Ng7
    """

    # Using the python-chess PGN parser is the most reliable method.
    pgn_io = io.StringIO(pgn_string)
    game = chess.pgn.read_game(pgn_io)

    # We navigate to the position just before black's 130th move.
    # White's 130th move is the 259th half-move in the game.
    node = game
    for _ in range(259):
        node = node.next()
    
    board = node.board()

    # The SAN move for Qe6 is just 'Qe6'
    blunder_move = board.parse_san('Qe6')
    
    # Candidate moves from the question
    # Maps letters to SAN moves
    candidate_moves = {
        'A': 'Qa1', 'B': 'Qa7', 'C': 'Qg2', 'D': 'Qf2',
        'E': 'Qb2', 'F': 'Qd2', 'G': 'Qa6', 'H': 'Qh2',
        'I': 'Qa8', 'J': 'Qa5', 'K': 'Qa4', 'L': 'Qc2',
        'M': 'Qe2', 'N': 'Qa3'
    }

    # Please make sure the path to the stockfish executable is correct for your system.
    # For Linux, it can often be installed via a package manager (e.g., `sudo apt-get install stockfish`)
    # and might be found at "/usr/games/stockfish". For Windows/macOS, you may need to download it
    # from the official Stockfish website and provide the full path to the executable.
    stockfish_path = "/usr/games/stockfish"
    try:
        engine = chess.engine.SimpleEngine.popen_uci(stockfish_path)
    except FileNotFoundError:
        print(f"Error: Stockfish engine not found at '{stockfish_path}'.")
        print("Please install Stockfish or update the 'stockfish_path' variable in the script.")
        return

    # Analyze the blunder
    board.push(blunder_move)
    info = engine.analyse(board, chess.engine.Limit(time=0.5))
    blunder_score = info["score"].white()
    print(f"Analysis of the actual move 130...Qe6:")
    print(f"Position is winning for White. Score: {blunder_score}\n")
    board.pop() # Revert move

    print("Analyzing alternative Queen moves...")
    best_move_info = {'choice': None, 'move': None, 'score': None}
    
    for choice, move_san in candidate_moves.items():
        try:
            move = board.parse_san(move_san)
            # Only analyze legal moves
            if move in board.legal_moves:
                board.push(move)
                info = engine.analyse(board, chess.engine.Limit(time=0.5))
                score = info["score"].white()
                print(f"  - Move {choice} ({move_san}): Score = {score}")
                
                # We are looking for the score closest to a draw (0.0)
                if best_move_info['score'] is None or abs(score.score(mate_score=10000)) < abs(best_move_info['score'].score(mate_score=10000)):
                    best_move_info = {'choice': choice, 'move': move_san, 'score': score}
                board.pop()

        except chess.InvalidMoveError:
             # This move is not legal from the current position
            pass
            
    engine.quit()

    if best_move_info['choice']:
        print(f"\nConclusion:")
        print(f"The move played, Qe6, led to a winning position for White (Score: {blunder_score}).")
        print(f"The best alternative is {best_move_info['move']}, which leads to a drawn position (Score: {best_move_info['score']}).")
        print(f"\nThe correct drawing move corresponds to answer choice {best_move_info['choice']}.")
        print(f"<<<{best_move_info['choice']}>>>")
    else:
        print("\nCould not find a valid drawing move among the choices.")

if __name__ == '__main__':
    find_drawing_move()