def is_valid(board, row, col, num):
    # Check if num is not in the current row
    if num in board[row]:
        return False
    
    # Check if num is not in the current column
    for r in range(5):
        if board[r][col] == num:
            return False
    
    return True

def solve_futoshiki(board):
    for row in range(5):
        for col in range(5):
            if board[row][col] == 0:
                for num in range(1, 6):
                    if is_valid(board, row, col, num):
                        board[row][col] = num
                        if solve_futoshiki(board):
                            return True
                        board[row][col] = 0
                return False
    return True

# Initial board setup
board = [
    [3, 2, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [5, 0, 3, 0, 0],
    [0, 5, 0, 2, 0],
    [0, 4, 1, 0, 0]
]

solve_futoshiki(board)

# Print the solved board
for row in board:
    print('   '.join(map(str, row)))