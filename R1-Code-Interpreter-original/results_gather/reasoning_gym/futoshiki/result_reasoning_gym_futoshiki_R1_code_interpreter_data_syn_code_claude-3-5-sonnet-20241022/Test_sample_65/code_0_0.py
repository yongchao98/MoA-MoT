def is_valid(board, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in board[row]:
        return False
    
    # Check column
    if num in [board[i][col] for i in range(5) if board[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    if col < 4:
        if h_constraints[row][col] == '<' and board[row][col+1] != 0:
            if not (num < board[row][col+1]):
                return False
        if h_constraints[row][col] == '>' and board[row][col+1] != 0:
            if not (num > board[row][col+1]):
                return False
    if col > 0:
        if h_constraints[row][col-1] == '<' and board[row][col-1] != 0:
            if not (board[row][col-1] < num):
                return False
        if h_constraints[row][col-1] == '>' and board[row][col-1] != 0:
            if not (board[row][col-1] > num):
                return False
    
    # Check vertical constraints
    if row < 4:
        if v_constraints[row][col] == 'v' and board[row+1][col] != 0:
            if not (num > board[row+1][col]):
                return False
        if v_constraints[row][col] == '^' and board[row+1][col] != 0:
            if not (num < board[row+1][col]):
                return False
    if row > 0:
        if v_constraints[row-1][col] == 'v' and board[row-1][col] != 0:
            if not (board[row-1][col] > num):
                return False
        if v_constraints[row-1][col] == '^' and board[row-1][col] != 0:
            if not (board[row-1][col] < num):
                return False
    
    return True

def find_empty(board):
    for i in range(5):
        for j in range(5):
            if board[i][j] == 0:
                return i, j
    return None

def solve_futoshiki(board, h_constraints, v_constraints):
    empty = find_empty(board)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 6):
        if is_valid(board, row, col, num, h_constraints, v_constraints):
            board[row][col] = num
            if solve_futoshiki(board, h_constraints, v_constraints):
                return True
            board[row][col] = 0
    
    return False

# Initialize the board
board = [
    [0, 0, 0, 0, 0],
    [0, 0, 4, 0, 0],
    [0, 0, 0, 3, 0],
    [0, 0, 3, 4, 0],
    [4, 0, 0, 0, 0]
]

# Horizontal constraints ('<' or '>' or ' ')
h_constraints = [
    [' ', '<', '>', '<', ' '],
    [' ', '<', '<', '>', ' '],
    [' ', ' ', '<', '<', ' '],
    [' ', '>', ' ', ' ', ' '],
    ['>', ' ', ' ', ' ', ' ']
]

# Vertical constraints ('^' or 'v' or ' ')
v_constraints = [
    ['^', ' ', 'v', ' ', ' '],
    [' ', ' ', ' ', ' ', '^'],
    [' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', 'v', ' '],
    [' ', ' ', ' ', ' ', ' ']
]

if solve_futoshiki(board, h_constraints, v_constraints):
    result = ""
    for i in range(5):
        row = ""
        for j in range(5):
            row += str(board[i][j])
            if j < 4:
                row += " " + h_constraints[i][j] + " "
        result += row + "\n"
        if i < 4:
            for j in range(5):
                if v_constraints[i][j] in ['^', 'v']:
                    result += v_constraints[i][j] + " " * 3
                else:
                    result += " " * 4
            result += "\n"
    print(result.rstrip())