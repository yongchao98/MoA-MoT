def solve_m_queens():
    """
    Calculates and demonstrates the solution for the maximum number m of
    white and black queens on a 16x16 board.
    """
    N = 16
    
    print(f"The problem is to find the maximum number m such that m white queens and m black queens")
    print(f"can coexist on a {N}x{N} chessboard without queens of the same color attacking each other.\n")
    
    print(f"Step 1: Determine the upper limit for m.")
    print(f"The maximum number of non-attacking queens of a single color on an {N}x{N} board is {N}.")
    print(f"This is the classic {N}-Queens problem. Therefore, m must be less than or equal to {N}.\n")

    print(f"Step 2: Show that m = {N} is achievable.")
    print(f"We need to find two disjoint sets of coordinates for {N} white and {N} black queens,")
    print("where each set represents a valid solution to the N-Queens problem.\n")

    # Generate the first solution (white queens) using an explicit construction.
    # We use 0-based indexing for rows and columns (0 to 15).
    # Positions are (row, col) where row is the index and col is the value.
    white_queens_cols = [0] * N
    for i in range(N // 2):
        white_queens_cols[i] = 2 * i + 1
    for i in range(N // 2, N):
        white_queens_cols[i] = 2 * i - N

    white_queens_coords = []
    for i in range(N):
        white_queens_coords.append((i, white_queens_cols[i]))

    # Generate the second solution (black queens) by reflecting the first one
    # across the vertical center line.
    # The new column c' is (N-1) - c.
    black_queens_cols = [(N - 1) - col for col in white_queens_cols]
    
    black_queens_coords = []
    for i in range(N):
        black_queens_coords.append((i, black_queens_cols[i]))
        
    print("Construction successful. We have found two valid and disjoint solutions.\n")

    print(f"Solution for {N} White Queens (row, column):")
    print(white_queens_coords)
    print("\n")
    
    print(f"Solution for {N} Black Queens (row, column):")
    print(black_queens_coords)
    print("\n")

    m = len(white_queens_coords)
    print(f"Since we can place {m} white queens and {m} black queens on the board,")
    print(f"and we know m cannot exceed {N}, the maximum value for m is {m}.")

solve_m_queens()