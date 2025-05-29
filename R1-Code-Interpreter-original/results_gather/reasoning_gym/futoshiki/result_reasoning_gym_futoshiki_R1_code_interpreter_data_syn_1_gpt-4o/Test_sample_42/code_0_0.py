def is_valid(board, row, col, num):
    # Check if num is not in the current row
    if num in board[row]:
        return False
    
    # Check if num is not in the current column
    for r in range(8):
        if board[r][col] == num:
            return False
    
    # Check inequality constraints
    # Vertical constraint between row 5 and 6
    if row == 4 and col == 4 and num <= board[5][4]:
        return False
    if row == 5 and col == 4 and num >= board[4][4]:
        return False
    
    return True

def solve_futoshiki(board):
    for row in range(8):
        for col in range(8):
            if board[row][col] == 0:  # Empty cell
                for num in range(1, 9):
                    if is_valid(board, row, col, num):
                        board[row][col] = num
                        if solve_futoshiki(board):
                            return True
                        board[row][col] = 0
                return False
    return True

# Initial board setup
board = [
    [8, 7, 4, 5, 0, 0, 1, 0],
    [0, 1, 2, 0, 0, 4, 0, 0],
    [2, 0, 0, 0, 0, 0, 6, 3],
    [1, 5, 0, 6, 0, 0, 3, 0],
    [5, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 8, 0, 5],
    [0, 8, 0, 0, 0, 0, 5, 4],
    [4, 0, 5, 7, 6, 0, 0, 2]
]

solve_futoshiki(board)

# Print the solved board
for row in board:
    print('   '.join(str(num) if num != 0 else '_' for num in row))