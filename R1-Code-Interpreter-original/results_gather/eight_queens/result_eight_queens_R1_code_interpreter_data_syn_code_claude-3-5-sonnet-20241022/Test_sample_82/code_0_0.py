def is_safe(board, row, col, forbidden):
    # Check if position is forbidden
    if (row, col) in forbidden:
        return False
    
    # Check row
    for j in range(8):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(8):
        if board[i][col] == 1:
            return False
    
    # Check diagonal (top-left to bottom-right)
    for i, j in zip(range(row-1, -1, -1), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row+1, 8), range(col+1, 8)):
        if board[i][j] == 1:
            return False
    
    # Check diagonal (top-right to bottom-left)
    for i, j in zip(range(row-1, -1, -1), range(col+1, 8)):
        if board[i][j] == 1:
            return False
    for i, j in zip(range(row+1, 8), range(col-1, -1, -1)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_queens(board, forbidden, queens_placed=0, start_row=0):
    if queens_placed == 8:
        return True, board
    
    for row in range(start_row, 8):
        for col in range(8):
            if board[row][col] == 0 and is_safe(board, row, col, forbidden):
                board[row][col] = 1
                success, result = solve_queens(board, forbidden, queens_placed + 1, row + 1)
                if success:
                    return True, result
                board[row][col] = 0
    
    return False, board

# Initialize board with given queen and forbidden positions
board = [[0]*8 for _ in range(8)]
board[0][2] = 1  # Initial queen
forbidden = {(3,1), (4,3), (7,6)}  # Forbidden positions

success, final_board = solve_queens(board, forbidden, 1)

if success:
    # Get queen positions
    positions = []
    for i in range(8):
        for j in range(8):
            if final_board[i][j] == 1:
                positions.append(f"{i} {j}")
    print(", ".join(positions))