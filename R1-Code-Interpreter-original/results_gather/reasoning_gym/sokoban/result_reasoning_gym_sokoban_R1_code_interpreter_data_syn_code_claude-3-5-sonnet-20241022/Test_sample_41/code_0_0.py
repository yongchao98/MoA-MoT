from collections import deque
import copy

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def is_goal_state(board):
    box_count = 0
    for row in board:
        for cell in row:
            if cell in ['$']:
                box_count += 1
    return box_count == 4

def get_valid_moves(board, player_pos):
    moves = []
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    for move, (dy, dx) in directions.items():
        new_y, new_x = player_pos[0] + dy, player_pos[1] + dx
        
        if board[new_y][new_x] == '+':
            continue
            
        if board[new_y][new_x] in ['@', '$']:
            box_y, box_x = new_y + dy, new_x + dx
            if board[box_y][box_x] not in ['+', '@', '$']:
                moves.append(move)
        else:
            moves.append(move)
    
    return moves

def make_move(board, move):
    new_board = [list(row) for row in board]
    player_pos = get_player_pos(new_board)
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    dy, dx = directions[move]
    
    py, px = player_pos
    new_py, new_px = py + dy, px + dx
    
    # Moving box
    if new_board[new_py][new_px] in ['@', '$']:
        box_y, box_x = new_py + dy, new_px + dx
        if new_board[new_py][new_px] == '@':
            if new_board[box_y][box_x] == 'X':
                new_board[box_y][box_x] = '$'
            else:
                new_board[box_y][box_x] = '@'
        else:  # box on goal
            if new_board[box_y][box_x] == 'X':
                new_board[box_y][box_x] = '$'
            else:
                new_board[box_y][box_x] = '@'
        new_board[new_py][new_px] = '*' if new_board[new_py][new_px] == '@' else '%'
    
    # Moving player
    if new_board[py][px] == '*':
        new_board[py][px] = '-'
    else:  # player on goal
        new_board[py][px] = 'X'
    
    if new_board[new_py][new_px] in ['-', 'X']:
        new_board[new_py][new_px] = '*' if new_board[new_py][new_px] == '-' else '%'
    
    return [''.join(row) for row in new_board]

def solve_sokoban(initial_board):
    queue = deque([(initial_board, "")])
    visited = set([''.join(initial_board)])
    
    while queue:
        current_board, path = queue.popleft()
        
        if is_goal_state(current_board):
            return path
        
        player_pos = get_player_pos(current_board)
        valid_moves = get_valid_moves(current_board, player_pos)
        
        for move in valid_moves:
            new_board = make_move(current_board, move)
            board_str = ''.join(new_board)
            
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
    "++++++++",
]

solution = solve_sokoban(initial_board)
print(solution)