def solve_chess_puzzle():
    """
    This function prints the step-by-step solution for the mate-in-2 chess puzzle.
    """
    print("This is a mate-in-2 puzzle. White to move.")
    print("The solution involves a key first move that forces checkmate on the next turn against any defense.")
    print("\n-------------------------")
    print("THE SOLUTION")
    print("-------------------------")

    print("\nStep 1: White's Key Move")
    print("1. Rh3")
    print("White moves the Rook from f3 to h3. This move prepares multiple checkmating threats.")

    print("\nStep 2: Analyzing Black's Replies and White's Mating Moves")

    print("\nVariation A: If Black captures the rook...")
    print("1. ...      Kxf3")
    print("The Black King takes the rook. White now delivers the final blow:")
    print("2. Ng6# (Checkmate)")
    print("The Knight from e7 moves to g6, checking the king. The king on f3 has no escape squares.")

    print("\nVariation B: If Black's king moves to d4...")
    print("1. ...      Kd4")
    print("The Black King escapes to d4. White delivers checkmate with the queen:")
    print("2. Qd7# (Checkmate)")
    print("The Queen moves to d7, checking the king. The king on d4 is trapped by the queen, the rook on h3, and the pawn on b4.")
    
    print("\nVariation C: If Black tries any other move...")
    print("For any other Black move (e.g., 1... Qe6), White proceeds with the original threat:")
    print("2. Qa1# (Checkmate)")
    print("The Queen moves to a1, delivering checkmate on the long diagonal.")

if __name__ == '__main__':
    solve_chess_puzzle()