import sys

def solve_chess_puzzle():
    """
    This function analyzes the given chess position and prints the solution.
    """
    # Introduction to the strategy
    print("The initial position is as follows:")
    print("Black: Rook(A8, F8), Knight(B8), Bishop(B7, F6), Queen(E7), King(G8), Pawns(A7, C7, D7, G7, H7, B6, E6)")
    print("White: Rook(A1, H1), Knight(E4, E5), Bishop(D3), Queen(H5), King(E1), Pawns(A2, B2, C2, D4, F2, G2, H2)")
    print("\nWhite is to move. The position is set up for a decisive attack on the Black king.")
    print("The winning combination involves a queen sacrifice to break through the king's defenses.")
    
    # The mating sequence
    print("\n--- The Mating Sequence ---")
    
    # Move 1
    print("\nMove: 1")
    print("White plays: Qxh7+")
    print("Explanation: White captures the pawn on h7 with the Queen. This is a sacrifice, but it's a check.")
    print("The Black g7-pawn cannot recapture the Queen because it is pinned to the king by the White Bishop on d3.")
    
    print("\nBlack plays: Kxh7")
    print("Explanation: Black's only legal move is for the King to capture the Queen on h7.")
    
    # Move 2
    print("\nMove: 2")
    print("White plays: Ng6#")
    print("Explanation: The White Knight from e5 moves to g6, delivering a devastating double check.")
    print("  - Check 1: The Knight on g6 now directly checks the Black King on h7.")
    print("  - Check 2: The Knight moving away reveals a discovered check from the White Rook on h1.")
    
    # Conclusion
    print("\n--- Result ---")
    print("This is checkmate. The Black King has no safe squares to move to, and a double check cannot be blocked or the checking pieces captured in one move.")
    print("\nWhite delivers checkmate in 2 moves.")

# Execute the function to print the solution.
solve_chess_puzzle()
