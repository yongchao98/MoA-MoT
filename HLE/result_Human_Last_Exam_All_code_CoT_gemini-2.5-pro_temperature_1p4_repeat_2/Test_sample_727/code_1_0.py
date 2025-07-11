import chess
import chess.engine
import asyncio
import io
import shutil
import sys

async def find_drawing_chess_move():
    """
    Analyzes a specific chess position from Carlsen-Nepomniachtchi 2021 Game 6
    to find the drawing move for Black among a list of choices.
    """
    # Use the PGN to reach the position right before black's 130th move.
    pgn_to_position = "1.d4 Nf6 2. Nf3 d5 3. g3 e6 4. Bg2 Be7 5. O-O O-O 6. b3 c5 7. dxc5 Bxc5 8. c4 dxc4 9. Qc2 Qe7 10. Nbd2 Nc6 11. Nxc4 b5 12. Nce5 Nb4 13. Qb2 Bb7 14. a3 Nc6 15. Nd3 Bb6 16. Bg5 Rfd8 17. Bxf6 gxf6 18. Rac1 Nd4 19. Nxd4 Bxd4 20. Qa2 Bxg2 21. Kxg2 Qb7+ 22. Kg1 Qe4 23. Qc2 a5 24. Rfd1 Kg7 25. Rd2 Rac8 26. Qxc8 Rxc8 27. Rxc8 Qd5 28. b4 a4 29. e3 Be5 30. h4 h5 31. Kh2 Bb2 32. Rc5 Qd6 33. Rd1 Bxa3 34. Rxb5 Qd7 35. Rc5 e5 36. Rc2 Qd5 37. Rdd2 Qb3 38. Ra2 e4 39. Nc5 Qxb4 40. Nxe4 Qb3 41. Rac2 Bf8 42. Nc5 Qb5 43. Nd3 a3 44. Nf4 Qa5 45. Ra2 Bb4 46. Rd3 Kh6 47. Rd1 Qa4 48. Rda1 Bd6 49. Kg1 Qb3 50. Ne2 Qd3 51. Nd4 Kh7 52. Kh2 Qe4 53. Rxa3 Qxh4+ 54. Kg1 Qe4 55. Ra4 Be5 56. Ne2 Qc2 57. R1a2 Qb3 58. Kg2 Qd5+ 59. f3 Qd1 60. f4 Bc7 61. Kf2 Bb6 62. Ra1 Qb3 63. Re4 Kg7 64. Re8 f5 65. Raa8 Qb4 66. Rac8 Ba5 67. Rc1 Bb6 68. Re5 Qb3 69. Re8 Qd5 70. Rcc8 Qh1 71. Rc1 Qd5 72. Rb1 Ba7 73. Re7 Bc5 74. Re5 Qd3 75. Rb7 Qc2 76. Rb5 Ba7 77. Ra5 Bb6 78. Rab5 Ba7 79. Rxf5 Qd3 80. Rxf7+ Kxf7 81. Rb7+ Kg6 82. Rxa7 Qd5 83. Ra6+ Kh7 84. Ra1 Kg6 85. Nd4 Qb7 86. Ra2 Qh1 87. Ra6+ Kf7 88. Nf3 Qb1 89. Rd6 Kg7 90. Rd5 Qa2+ 91. Rd2 Qb1 92. Re2 Qb6 93. Rc2 Qb1 94. Nd4 Qh1 95. Rc7+ Kf6 96. Rc6+ Kf7 97. Nf3 Qb1 98. Ng5+ Kg7 99. Ne6+ Kf7 100. Nd4 Qh1 101. Rc7+ Kf6 102. Nf3 Qb1 103. Rd7 Qb2+ 104. Rd2 Qb1 105. Ng1 Qb4 106. Rd1 Qb3 107. Rd6+ Kg7 108. Rd4 Qb2+ 109. Ne2 Qb1 110. e4 Qh1 111. Rd7+ Kg8 112. Rd4 Qh2+ 113. Ke3 h4 114. gxh4 Qh3+ 115. Kd2 Qxh4 116. Rd3 Kf8 117. Rf3 Qd8+ 118. Ke3 Qa5 119. Kf2 Qa7+ 120. Re3 Qd7 121. Ng3 Qd2+ 122. Kf3 Qd1+ 123. Re2 Qb3+ 124. Kg2 Qb7 125. Rd2 Qb3 126. Rd5 Ke7 127. Re5+ Kf7 128. Rf5+ Ke8 129. e5 Qa2+ 130. Kh3"

    game_io = io.StringIO(pgn_to_position)
    game = chess.pgn.read_game(game_io)
    board = game.board()
    for move in game.mainline_moves():
        board.push(move)

    candidate_moves = {
        "A": "Qa1", "B": "Qa7", "C": "Qg2", "D": "Qf2", "E": "Qb2",
        "F": "Qd2", "G": "Qa6", "H": "Qh2", "I": "Qa8", "J": "Qa5",
        "K": "Qa4", "L": "Qc2", "M": "Qe2", "N": "Qa3",
    }
    
    # Check if stockfish is available
    stockfish_path = shutil.which("stockfish")
    if not stockfish_path:
        print("Stockfish engine not found in your system's PATH.")
        print("Please install Stockfish from https://stockfishchess.org/download/ and ensure it's in your PATH.")
        return

    try:
        transport, engine = await chess.engine.popen_uci(stockfish_path)

        drawish_move = None
        min_abs_score = float('inf')

        print("Analyzing candidate moves (evaluation from White's perspective, 0.00 is a draw):\n")
        
        for letter, san_move in candidate_moves.items():
            try:
                temp_board = board.copy()
                move = temp_board.parse_san(san_move)
                temp_board.push(move)

                # Analyze with a reasonable depth for this endgame
                info = await engine.analyse(temp_board, chess.engine.Limit(depth=22))
                # Evaluation from the current player's perspective (White)
                score_obj = info["score"].white()
                
                # Convert PovScore to centipawns
                if score_obj.is_mate():
                    score_cp = 10000 * (1 if score_obj.mate() > 0 else -1)
                else:
                    score_cp = score_obj.score()
                
                print(f"Move {san_move} ({letter}): Evaluation = {score_cp/100:.2f} pawns")

                if abs(score_cp) < min_abs_score:
                    min_abs_score = abs(score_cp)
                    drawish_move = (letter, san_move)

            except chess.IllegalMoveError:
                print(f"Move {san_move} ({letter}): Is an illegal move.")
        
        print("\n--- Conclusion ---")
        if drawish_move:
            letter, move = drawish_move
            print(f"The move that leads to the most balanced (drawing) position is {move}.")
            print(f"This is answer choice: {letter}")
        else:
            print("Could not determine the best move.")

    finally:
        await engine.quit()

if __name__ == "__main__":
    # In some environments like Jupyter, a new event loop is needed.
    if sys.platform.startswith('win') and sys.version_info >= (3, 8):
        asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    
    asyncio.run(find_drawing_chess_move())

<<<D>>>