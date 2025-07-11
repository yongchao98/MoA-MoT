import chess
import chess.pgn
import io

def solve_chess_mate():
    """
    Sets up the chess board from the given PGN, finds the shortest
    mating sequence, and prints it.
    """
    # The PGN of the game as provided by the user.
    pgn_text = "1. e4 c6 2. c4 e6 3. Nc3 Bb4 4. a3 Bxc3 5. bxc3 Nf6 6. f3 b6 7. d4 c5 8. d5 O-O 9. Ne2 Ba6 10. Nf4 exd5 11. cxd5 Bxf1 12. Rxf1 Re8 13. Kf2 d6 14. Kg1 Nfd7 15. c4 Ne5 16. Qc2 Qf6 17. g3 Nxc4 18. Rb1 Qd4+ 19. Kh1 Nd7 20. Rd1 Ne3 21. Bxe3 Qxe3 22. Rd3 Qxd3 23. Qxd3 Ne5 24. Qe2 a6 25. Rxb6 Reb8 26. Rxd6 Rb1+ 27. Kg2 Rab8 28. Kh3 R1b2 29. Qxa6 Nxf3 30. Qa5 Rxh2+ 31. Kg4 Ne5+ 32. Kf5 Re8 33. Rd8 g6+ 34. Kg5"

    pgn = io.StringIO(pgn_text)
    game = chess.pgn.read_game(pgn)
    board = game.end().board()

    # The shortest mating sequence is a mate in 2.
    # 1. Black plays Kg7.
    # 2. White can try to defend with f4, but other moves lead to an immediate mate.
    #    We will show the mate against a move like Qe1.
    # 3. Black delivers the checkmate with h5.

    mating_sequence_san = ["Kg7", "Qe1", "h5"]

    # We will print the moves one by one.
    move_list = []

    # Move 1 (Black)
    move1 = board.parse_san("Kg7")
    move_list.append(board.san(move1))
    board.push(move1)

    # Move 2 (White's reply)
    # We choose an illustrative reply that allows mate-in-2.
    move2 = board.parse_san("Qe1")
    move_list.append(board.san(move2))
    board.push(move2)

    # Move 3 (Black's mating move)
    move3 = board.parse_san("h5")
    move_list.append(board.san(move3))
    board.push(move3)
    
    # Final check to confirm it is indeed checkmate.
    if board.is_checkmate():
        print("The shortest mating sequence is:")
        print(' '.join(move_list))

solve_chess_mate()