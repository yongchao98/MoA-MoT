def solve_chess_puzzle():
    """
    This function provides a step-by-step solution to the given chess problem,
    explaining the forced checkmate sequence for White.
    """
    # Introduction to the solution
    print("White can force a checkmate in 3 moves.")
    print("The winning combination involves a knight sacrifice to expose the Black king.")
    print("\n--- The Mating Sequence ---\n")

    # Detailed move-by-move breakdown
    print("Move 1 (White): Knight moves from E5 to F7, delivering a check.")
    print("   - Black's only legal move is to capture the attacking knight with the Rook on F8.")
    print("   - Full move: 1. Nf7+ Rxf7\n")

    print("Move 2 (White): Queen moves from H5 to H7, delivering another check.")
    print("   - Black's king is forced to move from G8 to F8.")
    print("   - Full move: 2. Qh7+ Kf8\n")

    print("Move 3 (White): Queen moves from H7 to H8, delivering checkmate.")
    print("   - The Black king on F8 is in check and has no escape squares.")
    print("   - The White Queen is protected by the Bishop on D3, so it cannot be captured.")
    print("   - Final move: 3. Qh8#\n")

    # Summary that fulfills the "equation" requirement
    print("Final Equation of Moves: (1. Nf7+) + (2. Qh7+) + (3. Qh8#)")
    print("\nTotal number of moves for White to checkmate: 3")

# Execute the function to display the solution
solve_chess_puzzle()
