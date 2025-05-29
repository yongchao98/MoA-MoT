def is_safe(board, row, col, queens):
    # Check row
    for j in range(8):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(8):
        if board[i][col] == 1:
            return False
    
    # Check diagonal (top-left to bottom-right)
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, 8), range(col, 8)):
        if board[i][j] == 1:
            return False
    
    # Check diagonal (top-right to bottom-left)
    for i, j in zip(range(row, -1, -1), range(col, 8)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row, 8), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_queens(board, queens_placed=0, queens=[]):
    if queens_placed == 8:
        return True, queens
    
    # Skip positions where queens are already placed
    start_row = 0
    while start_row < 8:
        found_queen = False
        for col in range(8):
            if board[start_row][col] == 1:
                found_queen = True
                break
        if not found_queen:
            break
        start_row += 1
    
    if start_row >= 8:
        return False, queens
    
    # Try placing queen in each column of the current row
    for col in range(8):
        if is_safe(board, start_row, col, queens):
            board[start_row][col] = 1
            queens.append((start_row, col))
            success, result = solve_queens(board, queens_placed + 1, queens)
            if success:
                return True, result
            board[start_row][col] = 0
            queens.pop()
    
    return False, queens

# Initialize the board with given queens
board = [
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,1,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],
    [1,0,0,0,0,0,0,0]
]

# Initialize queens list with existing queens
queens = [(1,5), (7,0)]

# Solve the puzzle
success, final_queens = solve_queens(board, 2, queens)

# Sort and format the output
if success:
    final_queens.sort(key=lambda x: x[0])
    result = ", ".join(f"{row} {col}" for row, col in final_queens)
    print(f"<<{result}>>")
else:
    print("No solution found")