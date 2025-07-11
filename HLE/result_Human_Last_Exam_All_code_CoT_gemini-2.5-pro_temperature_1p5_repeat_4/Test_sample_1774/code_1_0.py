import chess
import chess.engine
import sys

def solve_chess_puzzle():
    """
    Analyzes a chess position to find the fastest checkmate for White.

    This function sets up a specific chess board position from a FEN string,
    then uses the Stockfish chess engine to find a forced mate sequence.
    It prints the moves and the final answer.
    """
    # FEN representation of the board state described by the user.
    # Board Position: rn3rk1/pbppq1pp/1p1ppb2/4N2Q/3PN3/3B4/PPP2PPP/R3K2R w KQkq - 0 1
    # White: K(e1), Q(h5), R(a1,h1), B(d3), N(e4,e5), P(a2,b2,c2,d4,f2,g2,h2)
    # Black: k(g8), q(e7), r(a8,f8), b(b7,f6), n(b8), p(a7,b6,c7,d7,e6,g7,h7)
    fen_string = "rn3rk1/pbppq1pp/1p1ppb2/4N2Q/3PN3/3B4/PPP2PPP/R3K2R w KQkq - 0 1"
    board = chess.Board(fen_string)

    stockfish_path = "stockfish"
    if sys.platform == "win32":
        stockfish_path += ".exe"

    engine = None
    try:
        # The Stockfish engine executable must be in your system's PATH
        # or you must provide the full path to it.
        engine = chess.engine.SimpleEngine.popen_uci(stockfish_path)

        # Analyze the board to find the best line, setting a depth limit
        # to ensure the engine finds the short mate quickly.
        info = engine.analyse(board, chess.engine.Limit(depth=20))

        # The score will be Mate(x) if a forced mate is found
        score = info.get("score").white()
        mate_in = score.mate()

        if mate_in:
            print(f"White has a forced checkmate in {mate_in} moves.")
            print("The winning sequence is:")

            # Get the principal variation (the sequence of best moves)
            pv = info.get("pv")
            temp_board = board.copy()
            move_number = 1
            for i in range(0, len(pv), 2):
                white_move = pv[i]
                # Print White's move with move number
                print(f"{move_number}. {temp_board.san(white_move)}", end=" ")
                temp_board.push(white_move)

                # Check if there is a corresponding Black move
                if i + 1 < len(pv):
                    black_move = pv[i+1]
                    print(f"{temp_board.san(black_move)}", end=" ")
                    temp_board.push(black_move)
                
                print() # Newline for the next pair of moves
                move_number += 1
        else:
            print("No forced checkmate was found in the given position.")

    except (FileNotFoundError, chess.engine.EngineTerminatedError):
        print("Stockfish engine not found.")
        print("Please install Stockfish and ensure it's in your system's PATH.")
    finally:
        if engine:
            engine.quit()

if __name__ == '__main__':
    solve_chess_puzzle()
