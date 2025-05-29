from collections import deque

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def is_goal_state(board):
    goals = 0
    for row in board:
        goals += row.count('$')
    return goals == 4

def get_valid_moves(board, player_pos):
    moves = []
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    for move, (dy, dx) in directions.items():
        new_y, new_x = player_pos[0] + dy, player_pos[1] + dx
        
        # Check if new position is within bounds
        if new_y < 0 or new_y >= len(board) or new_x < 0 or new_x >= len(board[0]):
            continue
            
        if board[new_y][new_x] == '+':
            continue
            
        # Check if pushing a box
        if board[new_y][new_x] in ['@', '$']:
            box_y, box_x = new_y + dy, new_x + dx
            # Check if box position is within bounds
            if box_y < 0 or box_y >= len(board) or box_x < 0 or box_x >= len(board[0]):
                continue
            if board[box_y][box_x] not in ['+', '@', '$']:
                moves.append(move)
        else:
            moves.append(move)
    
    return moves

def make_move(board, move):
    board = [list(row) for row in board]
    player_pos = get_player_pos(board)
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    dy, dx = directions[move]
    
    py, px = player_pos
    new_py, new_px = py + dy, px + dx
    
    # Moving box
    if board[new_py][new_px] in ['@', '$']:
        box_y, box_x = new_py + dy, new_px + dx
        # Moving box from regular position
        if board[new_py][new_px] == '@':
            if board[box_y][box_x] == 'X':
                board[box_y][box_x] = '$'
            else:
                board[box_y][box_x] = '@'
        # Moving box from goal position
        else:
            if board[box_y][box_x] == 'X':
                board[box_y][box_x] = '$'
            else:
                board[box_y][box_x] = '@'
        board[new_py][new_px] = '*' if board[new_py][new_px] == '@' else '%'
    
    # Moving player
    if board[py][px] == '*':
        board[py][px] = '-'
    else:  # player on goal
        board[py][px] = 'X'
    
    if board[new_py][new_px] in ['-', 'X']:
        board[new_py][new_px] = '*' if board[new_py][new_px] == '-' else '%'
    
    return [''.join(row) for row in board]

def board_to_string(board):
    return '|'.join(board)

def solve_sokoban(initial_board):
    queue = deque([(initial_board, "")])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, path = queue.popleft()
        
        if is_goal_state(current_board):
            return path
        
        player_pos = get_player_pos(current_board)
        valid_moves = get_valid_moves(current_board, player_pos)
        
        for move in valid_moves:
            new_board = make_move(current_board, move)
            board_str = board_to_string(new_board)
            
            if board_str not in visited:
                visited.add(board_str)
                queue.append((new_board, path + move))
    
    return None

# Initial board
initial_board = [
    "++++++++",
    "+X-X+++",
    "+-X@-+",
    "+@@@X$-+",
    "+*----+",
    "++++++++"
]

solution = solve_sokoban(initial_board)
print(solution)