def find_player(board):
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def is_valid_move(board, x, y):
    return 0 <= x < len(board) and 0 <= y < len(board[0]) and board[x][y] != '+'

def is_complete(board):
    boxes = 0
    goals = 0
    boxes_on_goals = 0
    for row in board:
        for cell in row:
            if cell in ['@', '$']:
                boxes += 1
            if cell in ['X', '$', '%']:
                goals += 1
            if cell == '$':
                boxes_on_goals += 1
    return boxes == goals == boxes_on_goals

def try_move(board, move):
    board = [list(row) for row in board]
    player = find_player(board)
    if not player:
        return None

    px, py = player
    dx = -1 if move == 'U' else 1 if move == 'D' else 0
    dy = -1 if move == 'L' else 1 if move == 'R' else 0
    new_px, new_py = px + dx, py + dy

    if not is_valid_move(board, new_px, new_py):
        return None

    # Moving to empty space or goal
    if board[new_px][new_py] in ['-', 'X']:
        board[px][py] = '-' if board[px][py] == '*' else 'X'
        board[new_px][new_py] = '*' if board[new_px][new_py] == '-' else '%'
        return board

    # Pushing a box
    if board[new_px][new_py] in ['@', '$']:
        box_x, box_y = new_px + dx, new_py + dy
        if not is_valid_move(board, box_x, box_y) or board[box_x][box_y] in ['@', '$']:
            return None
        
        board[px][py] = '-' if board[px][py] == '*' else 'X'
        board[new_px][new_py] = '*' if board[new_px][new_py] == '@' else '%'
        board[box_x][box_y] = '@' if board[box_x][box_y] == '-' else '$'
        return board

    return None

def solve_sokoban():
    initial_board = [
        ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'],
        ['+', '-', '*', '-', '-', '-', '-', '-', 'X', '+'],
        ['+', '@', '@', '@', '@', '-', '-', '@', '-', '+'],
        ['+', '-', '-', '-', '-', 'X', '-', '$', '-', '+'],
        ['+', 'X', '-', '-', '-', '$', 'X', '-', '-', '+'],
        ['+', 'X', '-', '-', '-', '-', '-', '-', '-', '+'],
        ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
    ]
    
    # Modified solution after careful analysis
    solution = "RDDDDLLLUUURRDDDLLLUUURRDDDLLUUURRDDLLUURRDDLLUU"
    
    # Verify solution
    current_board = initial_board
    for move in solution:
        new_board = try_move(current_board, move)
        if not new_board:
            print("Invalid move:", move)
            return None
        current_board = new_board
    
    if is_complete(current_board):
        return solution
    return None

result = solve_sokoban()
print(result)