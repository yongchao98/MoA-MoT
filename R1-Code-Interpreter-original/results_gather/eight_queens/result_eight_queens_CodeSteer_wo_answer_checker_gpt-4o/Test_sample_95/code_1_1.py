def solve_n_queens_optimized(n):
    def backtrack(row, cols, diags1, diags2, positions):
        if row == n:
            return True
        
        for col in range(n):
            if not (cols & (1 << col)) and not (diags1 & (1 << (row - col + n - 1))) and not (diags2 & (1 << (row + col))):
                # Place the queen
                positions.append((row, col))
                # Update the constraints
                cols ^= (1 << col)
                diags1 ^= (1 << (row - col + n - 1))
                diags2 ^= (1 << (row + col))
                
                if backtrack(row + 1, cols, diags1, diags2, positions):
                    return True
                
                # Remove the queen and backtrack
                positions.pop()
                cols ^= (1 << col)
                diags1 ^= (1 << (row - col + n - 1))
                diags2 ^= (1 << (row + col))
        
        return False

    # Initialize constraints
    cols = 0
    diags1 = 0
    diags2 = 0
    positions = []
    
    # Create the board with constraints
    board = [
        [0, 'X', 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 'X', 0, 0],
        [0, 0, 0, 0, 0, 0, 'X', 0]
    ]
    
    # Pre-place queens based on the initial board
    for r in range(n):
        for c in range(n):
            if board[r][c] == 1:
                positions.append((r, c))
                cols |= (1 << c)
                diags1 |= (1 << (r - c + n - 1))
                diags2 |= (1 << (r + c))
    
    if backtrack(0, cols, diags1, diags2, positions):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in positions) + ">>>")
    else:
        print("No solution found")

solve_n_queens_optimized(8)