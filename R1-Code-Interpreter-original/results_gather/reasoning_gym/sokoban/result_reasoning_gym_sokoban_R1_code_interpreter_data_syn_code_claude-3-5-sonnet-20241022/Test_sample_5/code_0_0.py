from collections import deque
import copy

def is_win(board):
    # Check if all boxes are on goals
    goals = [(i, j) for i in range(len(board)) for j in range(len(board[0])) if board[i][j] in ['X', '$', '%']]
    boxes = [(i, j) for i in range(len(board)) for j in range(len(board[0])) if board[i][j] in ['@', '$']]
    return all(any(b == g for b in boxes) for g in goals)

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def is_valid_move(board, player_pos, direction):
    rows, cols = len(board), len(board[0])
    py, px = player_pos
    dy, dx = direction
    new_py, new_px = py + dy, px + dx
    
    # Check if new position is within bounds and not a wall
    if not (0 <= new_py < rows and 0 <= new_px < cols) or board[new_py][new_px] == '+':
        return False
        
    # If moving to empty space or goal
    if board[new_py][new_px] in ['-', 'X']:
        return True
        
    # If pushing a box, check if space behind box is valid
    if board[new_py][new_px] in ['@', '$']:
        box_py, box_px = new_py + dy, new_px + dx
        if not (0 <= box_py < rows and 0 <= box_px < cols):
            return False
        return board[box_py][box_px] in ['-', 'X']
        
    return False

def make_move(board, player_pos, direction):
    new_board = [list(row) for row in board]
    py, px = player_pos
    dy, dx = direction
    new_py, new_px = py + dy, px + dx
    
    # Update player position
    is_on_goal = board[py][px] == '%'
    new_board[py][px] = 'X' if is_on_goal else '-'
    
    # If pushing box
    if new_board[new_py][new_px] in ['@', '$']:
        box_py, box_px = new_py + dy, new_px + dx
        is_box_on_goal = new_board[new_py][new_px] == '$'
        new_board[new_py][new_px] = '%' if is_box_on_goal else '*'
        new_board[box_py][box_px] = '$' if new_board[box_py][box_px] == 'X' else '@'
    else:
        new_board[new_py][new_px] = '%' if new_board[new_py][new_px] == 'X' else '*'
    
    return [''.join(row) for row in new_board]

def solve_sokoban(initial_board):
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    queue = deque([(initial_board, "", get_player_pos(initial_board))])
    seen = {tuple(initial_board)}
    
    while queue:
        board, path, player_pos = queue.popleft()
        
        if is_win(board):
            return path
            
        for move, dir_vec in directions.items():
            if is_valid_move(board, player_pos, dir_vec):
                new_board = make_move(board, player_pos, dir_vec)
                board_tuple = tuple(new_board)
                
                if board_tuple not in seen:
                    seen.add(board_tuple)
                    queue.append((new_board, path + move, (player_pos[0] + dir_vec[0], player_pos[1] + dir_vec[1])))
    
    return None

# Initial board
initial_board = [
    "+++++++++++",
    "+-X----*-+",
    "+-@-@--@-+",
    "+----@X--+",
    "+-$X-X+$-+",
    "++++++++++",
]

solution = solve_sokoban(initial_board)
print(solution)