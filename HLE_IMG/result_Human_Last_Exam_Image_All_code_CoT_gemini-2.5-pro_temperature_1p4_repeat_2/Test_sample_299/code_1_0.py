def solve_chess_mate_in_2():
    """
    This function prints the step-by-step solution for the provided chess puzzle.
    White to move and checkmate in two moves.
    """
    print("The checkmate in 2 moves is found with a quiet first move that restricts the black king.")
    print("Here is the solution:")
    
    # White's first move
    # The notation includes the piece, starting square, and ending square for clarity.
    print("\nMove 1: White plays Rh5-h3")
    
    # Explanation of the two main variations after Black's reply.
    print("\nThis single move forces checkmate on the next turn, regardless of Black's reply.")
    print("There are two main possibilities for Black:")
    
    print("\n-------------------------------------------------------------")
    print("Variation 1: If Black plays 1... Qxf3")
    print("White delivers checkmate with: 2. Ne7-g6#")
    print("The knight on g6 checks the king. All escape squares are covered:")
    print("- d4 and f5 are covered by the Queen on a4.")
    print("- d3 is covered by the Bishop on g1.")
    print("- e3 is now covered by the Rook on h3.")
    print("-------------------------------------------------------------")

    print("\nVariation 2: If Black makes any other move (e.g., 1... Qd6)")
    print("White delivers checkmate with: 2. Rf3-e3#")
    print("The rook on e3 checks the king. The king cannot escape:")
    print("- d4 and f5 are covered by the Queen on a4.")
    print("- d3 is covered by the Bishop on g1.")
    print("- The king cannot capture the rook (Kxe3) because it would be in check from the Rook on h3.")
    print("-------------------------------------------------------------")

solve_chess_mate_in_2()