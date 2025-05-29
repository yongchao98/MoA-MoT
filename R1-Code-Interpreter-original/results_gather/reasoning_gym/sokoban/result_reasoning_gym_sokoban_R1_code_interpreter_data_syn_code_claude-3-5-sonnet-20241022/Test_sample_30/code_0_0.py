from collections import deque
from heapq import heappush, heappop
import copy

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def get_positions(board, char_list):
    positions = []
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in char_list:
                positions.append((i, j))
    return positions

def is_valid_move(board, player_pos, direction):
    rows, cols = len(board), len(board[0])
    new_pos = (player_pos[0] + direction[0], player_pos[1] + direction[1])
    
    if not (0 <= new_pos[0] < rows and 0 <= new_pos[1] < cols):
        return False
    
    if board[new_pos[0]][new_pos[1]] == '+':
        return False
        
    if board[new_pos[0]][new_pos[1]] in ['@', '$']:
        box_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        if not (0 <= box_pos[0] < rows and 0 <= box_pos[1] < cols):
            return False
        if board[box_pos[0]][box_pos[1]] in ['+', '@', '$']:
            return False
    return True

def make_move(board, player_pos, direction):
    new_board = [list(row) for row in board]
    new_pos = (player_pos[0] + direction[0], player_pos[1] + direction[1])
    
    # Handle player movement
    is_player_on_goal = board[player_pos[0]][player_pos[1]] == '%'
    new_board[player_pos[0]][player_pos[1]] = 'X' if is_player_on_goal else '-'
    
    # Handle box movement if necessary
    if board[new_pos[0]][new_pos[1]] in ['@', '$']:
        box_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        is_box_on_goal = board[new_pos[0]][new_pos[1]] == '$'
        new_board[new_pos[0]][new_pos[1]] = '%' if is_box_on_goal else '*'
        new_board[box_pos[0]][box_pos[1]] = '$' if board[box_pos[0]][box_pos[1]] == 'X' else '@'
    else:
        is_dest_goal = board[new_pos[0]][new_pos[1]] == 'X'
        new_board[new_pos[0]][new_pos[1]] = '%' if is_dest_goal else '*'
    
    return [''.join(row) for row in new_board], new_pos

def heuristic(board):
    boxes = get_positions(board, ['@', '$'])
    goals = get_positions(board, ['X', '$'])
    
    total_distance = 0
    for box in boxes:
        min_dist = min(manhattan_distance(box, goal) for goal in goals)
        total_distance += min_dist
    return total_distance

def solve_sokoban(initial_board):
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    initial_player = get_positions(initial_board, ['*', '%'])[0]
    
    queue = [(0, 0, initial_board, initial_player, [])]
    visited = set()
    
    while queue:
        _, cost, board, player_pos, path = heappop(queue)
        board_str = '\n'.join(board)
        
        if not get_positions(board, ['@']):  # All boxes are on goals
            return ''.join(path)
            
        if board_str in visited:
            continue
        visited.add(board_str)
        
        for move, dir_pos in directions.items():
            if is_valid_move(board, player_pos, dir_pos):
                new_board, new_player_pos = make_move(board, player_pos, dir_pos)
                new_cost = cost + 1
                priority = new_cost + heuristic(new_board)
                heappush(queue, (priority, new_cost, new_board, new_player_pos, path + [move]))
    
    return None

# Initial board
initial_board = [
    "++++++++",
    "+%@X---+",
    "+@-@@--+",
    "+------+",
    "+--XX--+",
    "+---@X-+",
    "++++++++",
]

solution = solve_sokoban(initial_board)
print(solution)