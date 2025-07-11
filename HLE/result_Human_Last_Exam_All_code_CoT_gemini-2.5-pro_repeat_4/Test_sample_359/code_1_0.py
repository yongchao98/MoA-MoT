# To run this script, you need to install the python-chess library.
# You can install it by running: pip install python-chess

import chess
import chess.variant

def solve_koth_puzzle():
    """
    This function analyzes the given King of the Hill chess puzzle,
    demonstrates the shortest winning line for White, and prints the result.
    """
    # The FEN string for the starting position of the puzzle.
    fen = "8/2k5/5pn1/1Pp1pNpp/3PP3/4K1B1/8/8 w - - 0 43"

    # Set up a King of the Hill board with the given FEN.
    board = chess.variant.KingOfTheHillBoard(fen)

    print("Initial Board State (FEN):")
    print(board.fen())
    print("-" * 20)

    # The optimal sequence of moves for a win in 4.
    # 1. dxc6+ Kb8 (Black avoids Kd8 which leads to mate in 2)
    # 2. Nxg6+ fxg6 (White removes the guardian of the center)
    # 3. Kf4 Ka7 (White's king moves towards e5, Black can't stop it)
    # 4. Ke5 (White's king reaches the center, winning the game)
    move_sequence_uci = ["d4c6", "c7b8", "f5g6", "f6g6", "e3f4", "b8a7", "f4e5"]
    
    # We will go through the moves one by one.
    # A full move in chess consists of a move by White and a move by Black.
    # The count starts at 1.
    move_count = 0
    # A ply is a half-move (one player's turn).
    for i, uci_move in enumerate(move_sequence_uci):
        ply = i + 1
        
        move = chess.Move.from_uci(uci_move)
        
        # Determine if it's a new move number
        if board.turn == chess.WHITE:
            move_count += 1
            move_san = board.san(move)
            print(f"Move {move_count}: White plays {move_san}")
        else:
            move_san = board.san(move)
            print(f"          Black plays {move_san}")

        board.push(move)

    print("-" * 20)
    print("Final Board State (FEN):")
    print(board.fen())
    print("\n")

    # Check the outcome of the game.
    # is_game_over() will be True due to the KotH win condition.
    if board.is_game_over():
        result = board.result()
        if result == "1-0":
            print(f"Result: White wins. The king is on e5, a central square.")
            print(f"The total number of moves for White to win is {move_count}.")
        else:
            print(f"The game ended, but not with a win for White. Result: {result}")
    else:
        print("The game is not over yet.")

if __name__ == '__main__':
    solve_koth_puzzle()
