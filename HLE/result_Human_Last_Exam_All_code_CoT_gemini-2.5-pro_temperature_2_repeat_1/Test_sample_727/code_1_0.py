import chess
import chess.engine
import sys
import os

def find_drawing_move():
    """
    Analyzes a chess position from Carlsen-Nepomniachtchi 2021, Game 6,
    to find the drawing move for Black among a list of choices.
    """
    
    # Step 1: Locate the Stockfish engine.
    # Common paths for Stockfish are checked.
    stockfish_paths = ["stockfish", "/usr/games/stockfish", "/usr/bin/stockfish"]
    engine_path = None
    for path in stockfish_paths:
        if os.path.exists(path) or any(os.access(os.path.join(p, path), os.X_OK) for p in os.environ["PATH"].split(os.pathsep)):
            engine_path = path
            break

    if not engine_path:
        print("Error: Stockfish engine not found.")
        print("Please install Stockfish from stockfishchess.org and ensure it's in your system's PATH.")
        sys.exit(1)

    try:
        # Step 2: Set up the board to the critical position.
        # The FEN string represents the board after White's move 130. Kh3.
        fen = "3k4/8/8/4PR2/8/6NK/q7/8 b - - 3 130"
        board = chess.Board(fen)
        
        # Step 3: Define candidate moves from the problem's answer choices.
        choices = {
            "A": "Qa1", "B": "Qa7", "C": "Qg2", "D": "Qf2", "E": "Qb2",
            "F": "Qd2", "G": "Qa6", "H": "Qh2", "I": "Qa8", "J": "Qa5",
            "K": "Qa4", "L": "Qc2", "M": "Qe2", "N": "Qa3"
        }
        
        # Add the actual losing move for comparison
        choices_with_blunder = choices.copy()
        choices_with_blunder["Blunder"] = "Qe6"

        results = {}
        print("Analyzing candidate moves for Black (evaluation from White's perspective)...")

        # Step 4: Analyze each move using the Stockfish engine.
        with chess.engine.SimpleEngine.popen_uci(engine_path) as engine:
            for letter, san_move in sorted(choices_with_blunder.items()):
                try:
                    move = board.parse_san(san_move)
                    temp_board = board.copy()
                    temp_board.push(move)
                    # Analyze the position for 0.5 seconds
                    info = engine.analyse(temp_board, chess.engine.Limit(time=0.5))
                    # Score is from the current player's (White's) perspective.
                    # A positive score means White is better. Black wants a score near 0.
                    score_cp = info["score"].white().score(mate_score=10000)
                    results[san_move] = score_cp
                    print(f"- Move {san_move:>4} ({letter}): Evaluation = {score_cp/100.0:+.2f}")
                except (chess.InvalidMoveError, chess.IllegalMoveError):
                    results[san_move] = "Illegal"
                    print(f"- Move {san_move:>4} ({letter}): Illegal Move")

        # Step 5: Identify the best drawing move.
        best_move_san = ""
        # A draw is ~0, so we look for the minimum absolute score.
        min_abs_score = float('inf')

        for san_move in choices.values():
            if results[san_move] != "Illegal":
                score = results[san_move]
                if abs(score) < min_abs_score:
                    min_abs_score = abs(score)
                    best_move_san = san_move
        
        print("\n--- Conclusion ---")
        losing_score = results.get("Qe6")
        if isinstance(losing_score, int):
            print(f"The blunder played in the game, Qe6, leads to an evaluation of +{losing_score/100.0:.2f}, which is a winning advantage for White.")
        
        best_score_cp = results.get(best_move_san, 0)
        best_score_pawn = best_score_cp / 100.0
        
        # Find the letter corresponding to the best move
        final_answer_letter = [letter for letter, san in choices.items() if san == best_move_san][0]

        print(f"\nThe best move to secure a draw is {best_move_san} (Choice {final_answer_letter}).")
        print("This move leads to an evaluation that is very close to zero.")
        print("\nThe final evaluation can be expressed as the following equation:")
        print(f"{best_score_cp} (centipawns) / 100 = {best_score_pawn:.2f} (pawns)")


    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    find_drawing_move()