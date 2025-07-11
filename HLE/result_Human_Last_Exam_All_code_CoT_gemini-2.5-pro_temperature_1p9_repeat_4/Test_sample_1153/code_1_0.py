# You need to install the stockfish library first.
# You can do this by running: pip install stockfish

from stockfish import Stockfish
import platform
import os

def find_stockfish_path():
    """Finds the path to the Stockfish executable."""
    # Check for the executable in common system paths
    for path in os.environ["PATH"].split(os.pathsep):
        # Look for 'stockfish' or 'stockfish.exe'
        for exe_name in ('stockfish', 'stockfish.exe'):
            exe_path = os.path.join(path, exe_name)
            if os.path.isfile(exe_path) and os.access(exe_path, os.X_OK):
                return exe_path
    
    # Pre-defined common paths for different OS
    system = platform.system()
    if system == "Windows":
        paths = ["C:/Stockfish/stockfish.exe"]
    elif system == "Linux":
        paths = ["/usr/games/stockfish", "/usr/bin/stockfish"]
    elif system == "Darwin": # macOS
        paths = ["/usr/local/bin/stockfish", "/opt/homebrew/bin/stockfish"]
    else:
        paths = []
        
    for path in paths:
        if os.path.exists(path):
            return path
            
    return None

# --- Main script ---

# FEN string for the chess position
fen_position = "8/P7/1np1p3/5k2/6p1/1P1NK3/8/8 w - - 0 1"

# Candidate moves from the user's question
# In Standard Algebraic Notation (SAN) and UCI notation
candidate_moves = {
    "A": {"san": "a8=Q", "uci": "a7a8q"},
    "B": {"san": "Nc5", "uci": "d3c5"},
    "C": {"san": "Kd4", "uci": "e3d4"},
    "D": {"san": "Kf2", "uci": "e3f2"},
    "E": {"san": "Nf4", "uci": "d3f4"},
    "F": {"san": "b4", "uci": "b3b4"},
}

stockfish_path = find_stockfish_path()

if not stockfish_path:
    print("Stockfish engine not found.")
    print("Please install Stockfish and add it to your system's PATH,")
    print("or specify the path manually in this script, for example:")
    print("stockfish = Stockfish(path=\"/path/to/stockfish\")")
else:
    try:
        # Initialize Stockfish with a reasonable thinking time (depth)
        stockfish = Stockfish(path=stockfish_path, depth=18)
        
        # Set the position
        stockfish.set_fen_position(fen_position)
        
        # Get the top 5 moves from Stockfish
        print("--- Stockfish Analysis ---")
        top_moves = stockfish.get_top_moves(5)

        if not top_moves:
             print("Could not get an evaluation from Stockfish.")
        else:
            best_move_uci = top_moves[0]['Move']
            best_move_eval_cp = top_moves[0]['Centipawn']
            best_move_eval_pawn = best_move_eval_cp / 100.0

            print(f"Position: {fen_position}\n")
            print("Stockfish's top moves:")
            for i, move_info in enumerate(top_moves):
                print(f"{i+1}. Move: {move_info['Move']}, Evaluation: {move_info['Centipawn']} centipawns")
            
            print("\n--- Conclusion ---")
            
            best_choice_letter = None
            best_choice_san = ""
            for letter, move_data in candidate_moves.items():
                if move_data['uci'] == best_move_uci:
                    best_choice_letter = letter
                    best_choice_san = move_data['san']
                    break
            
            if best_choice_letter:
                print(f"The best move is {best_choice_san}, which corresponds to option {best_choice_letter}.")
                print(f"The engine evaluates this move as a significant advantage for White.")
                print(f"Final Evaluation for {best_choice_san}: +{best_move_eval_pawn} pawn units")
            else:
                print(f"The best move according to Stockfish is {best_move_uci}, which is not in the options A-F.")
                print("Let's re-evaluate the best move among the given options.")
                best_option_move = ""
                best_option_eval = -99999
                # If the absolute best move isn't an option, check the options directly
                for letter, move_data in candidate_moves.items():
                    if stockfish.is_move_correct(move_data['uci']):
                        stockfish.make_moves_from_current_position([move_data['uci']])
                        evaluation = stockfish.get_evaluation()
                        print(f"Option {letter} ({move_data['san']}): Eval = {evaluation['value']} centipawns")
                        # We are white, so we want the max evaluation
                        if evaluation['value'] > best_option_eval:
                             best_option_eval = evaluation['value']
                             best_option_move = f"{letter} ({move_data['san']})"
                        stockfish.set_fen_position(fen_position) # reset position
                print(f"\nAmong the choices, the best one is {best_option_move}")
                

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure Stockfish is installed correctly and the path is accessible.")
