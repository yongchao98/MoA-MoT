def solve_chess_puzzle():
    """
    This function explains the solution to the chess puzzle and prints the mating sequence.
    """
    moves_to_mate = 3
    
    print(f"White can deliver checkmate in {moves_to_mate} moves.")
    print("The winning sequence involves a queen sacrifice followed by a decisive knight attack.")
    print("\nHere is the step-by-step mating sequence:")
    
    # Move 1
    print("\n---")
    print("Move 1 (White): Queen takes pawn on h7, check.")
    print("Black's only legal response is for the King to capture the Queen.")
    print("Equation: 1. Qxh7+ Kxh7")
    
    # Move 2
    print("\n---")
    print("Move 2 (White): Knight takes bishop on f6, creating a double check.")
    print("The King is checked by both the Knight on f6 and the Bishop on d3.")
    print("Black's King must move. We will follow the line where Black plays King to h6.")
    print("Equation: 2. Nxf6++ Kh6")

    # Move 3
    print("\n---")
    print("Move 3 (White): Knight from e5 moves to g4, delivering checkmate.")
    print("The Black King is trapped and cannot escape the attack.")
    print("Equation: 3. Nhg4#")
    print("---")

solve_chess_puzzle()
<<<3>>>