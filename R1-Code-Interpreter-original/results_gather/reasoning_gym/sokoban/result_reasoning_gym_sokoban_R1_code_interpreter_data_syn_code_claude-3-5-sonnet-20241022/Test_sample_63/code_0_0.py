from collections import deque
import copy

# Initial board
initial_board = [
    "+++++++++++",
    "+-$$X--@-+",
    "+-----X@-+",
    "+-@X@@-X-+",
    "+-@-----+",
    "+-@-X--$-+",
    "+-*X----+",
    "+-----X--+",
    "+---------+",
    "+++++++++++"
]

def is_box(char):
    return char in '@$'

def is_goal(char):
    return char in 'X$%'

def is_player(char):
    return char in '*%'

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if is_player(board[i][j]):
                return (i, j)
    return None

def get_board_string(board):
    return ''.join([''.join(row) for row in board])

def is_valid_move(board, player_pos, dx, dy):
    new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
    if board[new_x][new_y] == '+':
        return False
    if is_box(board[new_x][new_y]):
        box_x, box_y = new_x + dx, new_y + dy
        if board[box_x][box_y] in '+@$':
            return False
    return True

def make_move(board, direction):
    board = [list(row) for row in board]
    player_pos = get_player_pos(board)
    if not player_pos:
        return None
    
    dx, dy = 0, 0
    if direction == 'U': dx = -1
    elif direction == 'D': dx = 1
    elif direction == 'L': dy = -1
    elif direction == 'R': dy = 1
    
    if not is_valid_move(board, player_pos, dx, dy):
        return None
        
    px, py = player_pos
    new_px, new_py = px + dx, py + dy
    
    # Handle player movement
    player_on_goal = is_goal(board[px][py])
    new_pos_is_box = is_box(board[new_px][new_py])
    
    if new_pos_is_box:
        box_x, box_y = new_px + dx, new_py + dy
        box_to_goal = is_goal(board[box_x][box_y])
        board[box_x][box_y] = '$' if box_to_goal else '@'
        
    board[new_px][new_py] = '%' if is_goal(board[new_px][new_py]) else '*'
    board[px][py] = 'X' if player_on_goal else '-'
    
    return [''.join(row) for row in board]

def is_solved(board):
    for row in board:
        for char in row:
            if char in '@X':  # If there's a box not on goal or an empty goal
                return False
    return True

def solve_sokoban(initial_board):
    queue = deque([(initial_board, "")])
    visited = {get_board_string(initial_board)}
    
    while queue:
        current_board, path = queue.popleft()
        
        if is_solved(current_board):
            return path
            
        for move in 'UDLR':
            new_board = make_move(current_board, move)
            if new_board:
                board_string = get_board_string(new_board)
                if board_string not in visited:
                    visited.add(board_string)
                    queue.append((new_board, path + move))
    
    return None

solution = solve_sokoban(initial_board)
print(solution if solution else "No solution found")