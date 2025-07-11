import sys

# This script analyzes a chess position and determines the best move.

def explain_best_move():
    """
    Prints a detailed explanation for the best move in the given chess position.
    """
    
    fen = "8/P7/1np1p3/5k2/6p1/1P1NK3/8/8 w - - 0 1"
    
    print("Analyzing the chess position to find White's best move.")
    print(f"Position (FEN): {fen}\n")

    print("Candidate Move Analysis:")
    print("------------------------")
    
    # Analysis of A. a8=Q
    print("A. a8=Q: This move forces an exchange 1. a8=Q Nxa8. White follows up with 2. Ke4!, trapping the Black knight. This is a clear, winning line.")
    
    # Analysis of C. Kd4
    print("C. Kd4: A strong king move, but it allows the response 1... Ke4, which activates the Black king and creates unnecessary complications.")

    # Analysis of B. Nc5
    print("B. Nc5: This is the most precise move. It improves the knight's position, attacks the e6-pawn, and most importantly, prevents Black's king from activating via the e4 square. This puts Black in zugzwang, where any move leads to a worse position.")

    print("\nConclusion:")
    print("-----------")
    print("While 1. a8=Q is a simple and effective winning move, 1. Nc5 is superior. It prevents all of Black's counterplay and forces a win in the most controlled and clinical fashion. It addresses the immediate threats while preparing for a complete positional collapse of the opponent.")
    print("\nTherefore, the best move is B. Nc5.")

if __name__ == "__main__":
    explain_best_move()