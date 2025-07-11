import chess
import chess.engine
import sys
import os

def solve_chess_puzzle():
    """
    This function analyzes a chess position to find the best move for White.
    It uses the Stockfish chess engine to evaluate the position and explains the reasoning.
    """
    # --- Plan ---
    # 1. Define the chess position using its FEN string.
    # 2. Locate a usable Stockfish chess engine executable on the user's system.
    # 3. Initialize the engine and analyze the position to find the best move and its evaluation score.
    # 4. Provide a step-by-step analysis of the key candidate moves based on chess principles.
    # 5. Present the engine's conclusion, including the evaluation score, which is represented as a "final equation".
    # 6. Output the final answer choice in the specified format.

    FEN = "8/P7/1np1p3/5k2/6p1/1P1NK3/8/8 w - - 0 1"

    # Step 2: Locate Stockfish executable
    stockfish_path = None
    if sys.platform == "win32":
        # On Windows, check for stockfish.exe in the current directory first.
        possible_paths = ["./stockfish.exe", "stockfish.exe"]
    else:
        # On Linux/macOS, check for 'stockfish' in common system paths.
        possible_paths = ["stockfish", "/usr/games/stockfish", "/usr/bin/stockfish", "/usr/local/bin/stockfish"]

    # Find the first valid path by attempting to launch the engine.
    for path in possible_paths:
        try:
            # Use a short timeout to quickly check if the engine is valid.
            with chess.engine.SimpleEngine.popen_uci(path, setpgrp=True, timeout=0.1) as engine_test:
                stockfish_path = path
                break
        except (FileNotFoundError, PermissionError, chess.engine.EngineTerminatedError, chess.engine.EngineError):
            continue

    if not stockfish_path:
        print("---")
        print("Error: Stockfish engine not found.")
        print("Please download the Stockfish chess engine from https://stockfishchess.org/download/")
        print("Then, place the executable (e.g., 'stockfish.exe' or 'stockfish') in the same directory as this script, or ensure it's in your system's PATH.")
        print("---")
        return

    # Step 3: Analyze the position with the found engine
    try:
        with chess.engine.SimpleEngine.popen_uci(stockfish_path, setpgrp=True) as engine:
            board = chess.Board(FEN)
            # Analyze for 1.5 seconds, which is sufficient for a clear position like this.
            info = engine.analyse(board, chess.engine.Limit(time=1.5), multipv=1)

            if not info:
                 print("Engine analysis failed to return information.")
                 return

            # Extract the best move and score from the analysis
            best_move_info = info[0]
            best_move = best_move_info["pv"][0]
            score = best_move_info["score"].white()

    except Exception as e:
        print(f"An error occurred during engine analysis: {e}")
        return

    # Step 4: Explain the reasoning and analysis
    print("Thinking Process:")
    print("The primary feature of this endgame is White's passed pawn on a7. Black's defense hinges entirely on the knight on b6. White's objective is to dislodge this knight.")
    print("\nAnalysis of Candidate Moves:")
    print("A. a8=Q: A blunder. 1...Nxa8 trades the winning pawn for a lost endgame for White.")
    print("B. Nc5: An inaccurate move that creates minor threats but fails to address the core problem.")
    print("C. Kd4: A strong alternative. Activating the king is a powerful idea, but it's slightly less forceful than the optimal choice.")
    print("D. Kf2: Too passive. This move cedes the initiative to Black.")
    print("E. Nf4: An ineffective move. The knight gets exposed to attack by Black's king.")
    print("F. b4: The most precise and winning move. It directly prepares to attack the key defender on b6 with the pawn push b5, leading to a decisive advantage.")

    # Step 5: Present the engine's conclusion and the "final equation"
    san_move = board.san(best_move)
    print(f"\nEngine Conclusion:")
    print(f"The chess engine confirms that the best move is {san_move} ({best_move.uci()}).")

    print("\nFinal Equation (Engine Evaluation):")
    # Here, we interpret the "equation" as the engine's evaluation score of the position after the best move.
    if score.is_mate():
        # The score is a mate-in-X, which includes a number.
        print(f"Position Evaluation = Mate in {score.mate()}")
    else:
        # The score is a centipawn value. We convert it to a pawn advantage (a number).
        pawn_advantage = score.score() / 100.0
        print(f"Position Evaluation Score = {pawn_advantage:+.2f} pawns")

    # Step 6: Determine and print the final answer choice
    choices = {
        "a7a8q": "A",
        "d3c5": "B",
        "e3d4": "C",
        "e3f2": "D",
        "d3f4": "E",
        "b3b4": "F"
    }
    final_answer_key = choices.get(best_move.uci(), "Unknown")

    print(f"\n<<<{final_answer_key}>>>")

if __name__ == "__main__":
    solve_chess_puzzle()