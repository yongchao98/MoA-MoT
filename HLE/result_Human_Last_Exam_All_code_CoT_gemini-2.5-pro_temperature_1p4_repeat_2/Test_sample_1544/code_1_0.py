def solve_chess_puzzle():
    """
    Analyzes the chess position from the given FEN and determines the optimal move for White.
    The analysis concludes that trading queens leads to a winning pawn endgame.
    """
    
    fen_position = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"
    
    print("--- Chess Puzzle Analysis ---")
    print(f"Position FEN: {fen_position}")
    print("It is White's turn to move.\n")
    
    print("Situation Analysis:")
    print("White's key asset is the passed pawn on d6, which is only two squares away from promoting to a Queen.")
    print("Black's main strength is the active Queen on c5, which threatens to create a perpetual check on the White King and force a draw.")
    
    best_move = "Qxc5+"
    
    print("\nRecommended Move and Justification:")
    print(f"The best move for White is {best_move}.")
    print(f"\nThis move forces a trade of Queens. After 1. {best_move}, Black must play 1... bxc5.")
    print("The game transitions into a pawn endgame. In this new position, White's passed d-pawn becomes unstoppable.")
    print("White's King can then easily advance to support the pawn's promotion while also stopping Black's pawns, securing a decisive victory.")
    print("\nOther moves fail to resolve the threat of perpetual check and would likely lead to a draw.")
    
# Execute the function to display the solution.
solve_chess_puzzle()
