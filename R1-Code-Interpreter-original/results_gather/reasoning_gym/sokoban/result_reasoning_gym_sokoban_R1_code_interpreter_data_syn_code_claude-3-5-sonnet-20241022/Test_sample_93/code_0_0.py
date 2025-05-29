from collections import deque
import copy

def get_initial_state():
    board = [
        ['+', '+', '+', '+', '+', '+', '+'],
        ['+', '+', '$', '-', '-', '-', '+'],
        ['+', '-', '-', '-', '@', '-', '+'],
        ['+', '-', 'X', '-', '-', '-', '+'],
        ['+', '@', '@', '-', '-', '-', '+'],
        ['+', '%', '$', 'X', '+', '+', '+'],
        ['+', '+', '+', '+', '+', '+', '+']
    ]
    return board

def get_player_positions(board):
    players = []
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                players.append((i, j))
    return players

def get_box_positions(board):
    boxes = []
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['@', '$']:
                boxes.append((i, j))
    return boxes

def is_goal_position(board, pos):
    return board[pos[0]][pos[1]] in ['X', '$', '%']

def board_to_string(board):
    return ''.join([''.join(row) for row in board])

def is_valid_move(board, player_pos, direction):
    new_pos = (player_pos[0] + direction[0], player_pos[1] + direction[1])
    if board[new_pos[0]][new_pos[1]] == '+':
        return False
    if board[new_pos[0]][new_pos[1]] in ['@', '$']:
        box_new_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        if board[box_new_pos[0]][box_new_pos[1]] in ['+', '@', '$']:
            return False
    return True

def make_move(board, player_pos, direction):
    new_board = copy.deepcopy(board)
    new_pos = (player_pos[0] + direction[0], player_pos[1] + direction[1])
    
    # Update player position
    if board[player_pos[0]][player_pos[1]] == '%':
        new_board[player_pos[0]][player_pos[1]] = 'X'
    else:
        new_board[player_pos[0]][player_pos[1]] = '-'
    
    if board[new_pos[0]][new_pos[1]] in ['@', '$']:
        box_new_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        if board[new_pos[0]][new_pos[1]] == '$':
            new_board[new_pos[0]][new_pos[1]] = 'X'
        else:
            new_board[new_pos[0]][new_pos[1]] = '-'
            
        if is_goal_position(board, box_new_pos):
            new_board[box_new_pos[0]][box_new_pos[1]] = '$'
        else:
            new_board[box_new_pos[0]][box_new_pos[1]] = '@'
    
    if is_goal_position(board, new_pos):
        new_board[new_pos[0]][new_pos[1]] = '%'
    else:
        new_board[new_pos[0]][new_pos[1]] = '*'
    
    return new_board

def is_solved(board):
    goals = 0
    boxes_on_goals = 0
    for row in board:
        for cell in row:
            if cell in ['X', '$', '%']:
                goals += 1
            if cell in ['$']:
                boxes_on_goals += 1
    return goals == boxes_on_goals

def solve_sokoban():
    initial_board = get_initial_state()
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    visited = set()
    queue = deque([(initial_board, "", get_player_positions(initial_board)[0])])
    
    while queue:
        current_board, moves, player_pos = queue.popleft()
        if is_solved(current_board):
            return moves
        
        board_state = board_to_string(current_board)
        if board_state in visited:
            continue
        visited.add(board_state)
        
        for move, dir_coords in directions.items():
            if is_valid_move(current_board, player_pos, dir_coords):
                new_board = make_move(current_board, player_pos, dir_coords)
                new_player_pos = (player_pos[0] + dir_coords[0], player_pos[1] + dir_coords[1])
                queue.append((new_board, moves + move, new_player_pos))
    
    return None

solution = solve_sokoban()
print(solution)