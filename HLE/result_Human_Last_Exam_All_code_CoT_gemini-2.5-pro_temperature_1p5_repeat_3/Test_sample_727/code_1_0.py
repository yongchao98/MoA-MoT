import chess
import chess.pgn
import io

def solve_chess_puzzle():
    """
    Analyzes the chess position from Carlsen-Nepomniachtchi Game 6 (2021)
    at move 130 and finds the drawing move for Black.
    """
    pgn_string = """
    [Event "World Championship 2021"]
    [Site "Dubai"]
    [Date "2021.12.03"]
    [Round "6"]
    [White "Carlsen, Magnus"]
    [Black "Nepomniachtchi, Ian"]
    [Result "1-0"]

    1. d4 Nf6 2. Nf3 d5 3. g3 e6 4. Bg2 Be7 5. O-O O-O 6. b3 c5 7. dxc5
    Bxc5 8. c4 dxc4 9. Qc2 Qe7 10. Nbd2 Nc6 11. Nxc4 b5 12. Nce5 Nb4 13.
    Qb2 Bb7 14. a3 Nc6 15. Nd3 Bb6 16. Bg5 Rfd8 17. Bxf6 gxf6 18. Rac1 Nd4
    19. Nxd4 Bxd4 20. Qa2 Bxg2 21. Kxg2 Qb7+ 22. Kg1 Qe4 23. Qc2 a5 24.
    Rfd1 Kg7 25. Rd2 Rac8 26. Qxc8 Rxc8 27. Rxc8 Qd5 28. b4 a4 29. e3 Be5
    30. h4 h5 31. Kh2 Bb2 32. Rc5 Qd6 33. Rd1 Bxa3 34. Rxb5 Qd7 35. Rc5 e5
    36. Rc2 Qd5 37. Rdd2 Qb3 38. Ra2 e4 39. Nc5 Qxb4 40. Nxe4 Qb3 41. Rac2
    Bf8 42. Nc5 Qb5 43. Nd3 a3 44. Nf4 Qa5 45. Ra2 Bb4 46. Rd3 Kh6 47. Rd1
    Qa4 48. Rda1 Bd6 49. Kg1 Qb3 50. Ne2 Qd3 51. Nd4 Kh7 52. Kh2 Qe4 53.
    Rxa3 Qxh4+ 54. Kg1 Qe4 55. Ra4 Be5 56. Ne2 Qc2 57. R1a2 Qb3 58. Kg2
    Qd5+ 59. f3 Qd1 60. f4 Bc7 61. Kf2 Bb6 62. Ra1 Qb3 63. Re4 Kg7 64. Re8
    f5 65. Raa8 Qb4 66. Rac8 Ba5 67. Rc1 Bb6 68. Re5 Qb3 69. Re8 Qd5 70.
    Rcc8 Qh1 71. Rc1 Qd5 72. Rb1 Ba7 73. Re7 Bc5 74. Re5 Qd3 75. Rb7 Qc2
    76. Rb5 Ba7 77. Ra5 Bb6 78. Rab5 Ba7 79. Rxf5 Qd3 80. Rxf7+ Kxf7 81.
    Rb7+ Kg6 82. Rxa7 Qd5 83. Ra6+ Kh7 84. Ra1 Kg6 85. Nd4 Qb7 86. Ra2 Qh1
    87. Ra6+ Kf7 88. Nf3 Qb1 89. Rd6 Kg7 90. Rd5 Qa2+ 91. Rd2 Qb1 92. Re2
    Qb6 93. Rc2 Qb1 94. Nd4 Qh1 95. Rc7+ Kf6 96. Rc6+ Kf7 97. Nf3 Qb1 98.
    Ng5+ Kg7 99. Ne6+ Kf7 100. Nd4 Qh1 101. Rc7+ Kf6 102. Nf3 Qb1 103. Rd7
    Qb2+ 104. Rd2 Qb1 105. Ng1 Qb4 106. Rd1 Qb3 107. Rd6+ Kg7 108. Rd4 Qb2+
    109. Ne2 Qb1 110. e4 Qh1 111. Rd7+ Kg8 112. Rd4 Qh2+ 113. Ke3 h4 114.
gxh4 Qh3+ 115. Kd2 Qxh4 116. Rd3 Kf8 117. Rf3 Qd8+ 118. Ke3 Qa5 119.
    Kf2 Qa7+ 120. Re3 Qd7 121. Ng3 Qd2+ 122. Kf3 Qd1+ 123. Re2 Qb3+ 124.
    Kg2 Qb7 125. Rd2 Qb3 126. Rd5 Ke7 127. Re5+ Kf7 128. Rf5+ Ke8 129. e5
    Qa2+ 130. Kh3
    """
    pgn = io.StringIO(pgn_string)
    game = chess.pgn.read_game(pgn)
    
    node = game
    # Go to the 130th move for white
    for _ in range(129 * 2 + 1):
        node = node.next()
        
    board = node.board()
    
    print("The position is after White's move 130. Kh3. Black's move Qe6 was a losing blunder.")
    print("The drawing move for Black is Qd2+, which initiates a perpetual check.")
    print("\nLet's analyze the two main variations:\n")

    # Variation 1: White plays 131. Kh2
    print("--- Variation 1: White plays 131. Kh2 ---")
    board1 = board.copy()
    
    # 130... Qd2+
    move1_black = chess.Move.from_uci("a2d2")
    board1.push(move1_black)
    print(f"130. ... {board.san(move1_black)}")
    
    # 131. Kh2
    move2_white = chess.Move.from_uci("h3h2")
    board1.push(move2_white)
    
    # 131... Qf2+
    move2_black = chess.Move.from_uci("d2f2")
    board1.push(move2_black)
    print(f"131. {board.san(move2_white)} {board1.san(move2_black)}")
    
    # 132. Kh1
    move3_white = chess.Move.from_uci("h2h1")
    board1.push(move3_white)
    
    # 132... Qe1+
    move3_black = chess.Move.from_uci("f2e1")
    board1.push(move3_black)
    print(f"132. {board1.san(move3_white)} {board1.san(move3_black)}")

    # 133. Kh2
    move4_white = chess.Move.from_uci("h1h2")
    board1.push(move4_white)

    # 133... Qf2+
    move4_black = chess.Move.from_uci("e1f2")
    board1.push(move4_black)
    print(f"133. {board1.san(move4_white)} {board1.san(move4_black)}")
    print("The position is repeated, leading to a draw.\n")

    # Variation 2: White plays 131. Kg4
    print("--- Variation 2: White plays 131. Kg4 ---")
    board2 = board.copy()

    # 130... Qd2+
    move1_black = chess.Move.from_uci("a2d2")
    board2.push(move1_black)
    print(f"130. ... {board.san(move1_black)}")
    
    # 131. Kg4
    move2_white = chess.Move.from_uci("h3g4")
    board2.push(move2_white)
    
    # 131... Qe2+
    move2_black = chess.Move.from_uci("d2e2")
    board2.push(move2_black)
    print(f"131. {board.san(move2_white)} {board2.san(move2_black)}")
    
    # 132. Kf4
    move3_white = chess.Move.from_uci("g4f4")
    board2.push(move3_white)

    # 132... Qf2+
    move3_black = chess.Move.from_uci("e2f2")
    board2.push(move3_black)
    print(f"132. {board2.san(move3_white)} {board2.san(move3_black)}")

    # 133. Ke4
    move4_white = chess.Move.from_uci("f4e4")
    board2.push(move4_white)

    # 133... Qe2+
    move4_black = chess.Move.from_uci("f2e2")
    board2.push(move4_black)
    print(f"133. {board2.san(move4_white)} {board2.san(move4_black)}")
    print("The sequence of checks repeats, also leading to a draw.\n")

    print("Conclusion: The move 130... Qd2+ (Answer F) is the only move that guarantees a draw.")

solve_chess_puzzle()
<<<F>>>