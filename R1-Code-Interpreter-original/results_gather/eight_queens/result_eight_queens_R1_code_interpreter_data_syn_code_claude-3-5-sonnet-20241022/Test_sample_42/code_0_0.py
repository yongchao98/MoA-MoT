def is_safe(board, row, col, n):
    # Check row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check diagonal (top-left to bottom-right)
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, n, 1), range(col, n, 1)):
        if board[i][j] == 1:
            return False
    
    # Check diagonal (top-right to bottom-left)
    for i, j in zip(range(row, -1, -1), range(col, n, 1)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, n, 1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens(board, n):
    queens_positions = []
    # First, get existing queens positions
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queens_positions.append((i, j))
    
    def solve(row, queens_left):
        if queens_left == 0:
            return True
        
        if row >= n:
            return False
        
        # Skip rows that already have queens
        while any(pos[0] == row for pos in queens_positions):
            row += 1
            if row >= n:
                return False
        
        for col in range(n):
            if is_safe(board, row, col, n):
                board[row][col] = 1
                queens_positions.append((row, col))
                
                if solve(row + 1, queens_left - 1):
                    return True
                
                board[row][col] = 0
                queens_positions.pop()
        
        return False

    # Initialize board with given queens
    initial_board = [
        [0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0]
    ]
    
    if solve(0, 5):  # We need to place 5 more queens
        # Sort positions for consistent output
        queens_positions.sort()
        result = ", ".join(f"{pos[0]} {pos[1]}" for pos in queens_positions)
        print(f"<<{result}>>")
    else:
        print("No solution exists")

# Solve the puzzle
solve_n_queens([[0]*8 for _ in range(8)], 8)