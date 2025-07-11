import chess
import chess.variant

def solve_king_of_the_hill():
    """
    Analyzes a specific King of the Hill chess position to find how
    quickly White can force a win.
    """
    fen = "8/2k5/5pn1/1Pp1pNpp/3PP3/4K1B1/8/8 w - - 0 43"
    
    # In King of the Hill, a win can be achieved by moving one's king
    # to the center squares (d4, e4, d5, e5).
    #
    # Our analysis suggests a 2-move win for White:
    # 1. White plays dxc5, clearing the d4 square.
    # 2. Black makes a move.
    # 3. White plays Kd4, moving the king to the winning square.
    #
    # This script will verify that for EVERY possible Black response to 1. dxc5,
    # White can indeed win with 2. Kd4.

    board = chess.variant.KingOfTheHillBoard(fen)
    
    # Define White's first move (d4xc5) and second move (Ke3-d4)
    white_move_1 = chess.Move.from_uci("d4c5")
    white_move_2 = chess.Move.from_uci("e3d4")
    
    # Make White's first move
    board.push(white_move_1)
    
    # Assume the line is a forced win until a counterexample is found
    is_forced_win = True
    
    # Iterate through all of Black's legal replies
    black_replies = list(board.legal_moves)
    
    if not black_replies and not board.is_game_over():
        # This case is unlikely, but if black has no moves and it's not
        # checkmate/stalemate, something is wrong.
        is_forced_win = False

    for black_move in black_replies:
        temp_board = board.copy()
        temp_board.push(black_move)
        
        # After black's reply, check if white's winning move is legal and wins
        if white_move_2 in temp_board.legal_moves:
            temp_board.push(white_move_2)
            # is_variant_win() confirms a King of the Hill victory
            if not temp_board.is_variant_win():
                is_forced_win = False
                break
        else:
            # If the winning move is not legal after a certain Black reply,
            # then it is not a forced win.
            is_forced_win = False
            break

    # If the win is confirmed after checking all of black's replies, print the answer.
    if is_forced_win:
        # The number of moves White takes to win is 2.
        print(2)
    else:
        # This part should not be reached for the given FEN
        print("The hypothesized 2-move win could not be verified.")

solve_king_of_the_hill()