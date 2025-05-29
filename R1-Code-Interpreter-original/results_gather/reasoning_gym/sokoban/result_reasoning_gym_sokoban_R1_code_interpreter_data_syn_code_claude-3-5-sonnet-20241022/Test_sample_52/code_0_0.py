from collections import deque
import heapq
from copy import deepcopy

def parse_board(board_str):
    return [list(row.strip()) for row in board_str.strip().split('\n')]

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['@', '*']:
                return (i, j)
    return None

def get_boxes(board):
    boxes = set()
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['$', '@']:
                boxes.add((i, j))
    return boxes

def get_goals(board):
    goals = set()
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['X', '*']:
                goals.add((i, j))
    return goals

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def heuristic(boxes, goals):
    if not boxes or not goals:
        return 0
    total = 0
    boxes = list(boxes)
    goals = list(goals)
    for box in boxes:
        min_dist = float('inf')
        for goal in goals:
            dist = manhattan_distance(box, goal)
            min_dist = min(min_dist, dist)
        total += min_dist
    return total

def is_valid_move(board, player_pos, direction):
    new_pos = (player_pos[0] + direction[0], player_pos[1] + direction[1])
    if board[new_pos[0]][new_pos[1]] == '+':
        return False
    if board[new_pos[0]][new_pos[1]] in ['$', '@']:
        box_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        if board[box_pos[0]][box_pos[1]] in ['+', '$', '@']:
            return False
    return True

def make_move(board, player_pos, direction):
    new_board = [row[:] for row in board]
    new_pos = (player_pos[0] + direction[0], player_pos[1] + direction[1])
    
    # Update player position
    if board[player_pos[0]][player_pos[1]] == '*':
        new_board[player_pos[0]][player_pos[1]] = 'X'
    else:
        new_board[player_pos[0]][player_pos[1]] = '-'
    
    if board[new_pos[0]][new_pos[1]] == 'X':
        new_board[new_pos[0]][new_pos[1]] = '*'
    else:
        new_board[new_pos[0]][new_pos[1]] = '@'
    
    # If pushing a box
    if board[new_pos[0]][new_pos[1]] in ['$', '@']:
        box_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        if board[box_pos[0]][box_pos[1]] == 'X':
            new_board[box_pos[0]][box_pos[1]] = '$'
        else:
            new_board[box_pos[0]][box_pos[1]] = '@'
    
    return new_board

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_sokoban(initial_board):
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    initial_state = initial_board
    goals = get_goals(initial_state)
    visited = set()
    pq = [(0, '', initial_state)]
    
    while pq:
        _, path, current_state = heapq.heappop(pq)
        state_str = board_to_string(current_state)
        
        if state_str in visited:
            continue
        visited.add(state_str)
        
        # Check if won
        boxes = get_boxes(current_state)
        if boxes.issubset(goals):
            return path
        
        player_pos = get_player_pos(current_state)
        
        for move, dir_offset in directions.items():
            if is_valid_move(current_state, player_pos, dir_offset):
                new_state = make_move(current_state, player_pos, dir_offset)
                new_boxes = get_boxes(new_state)
                priority = len(path) + 1 + heuristic(new_boxes, goals)
                heapq.heappush(pq, (priority, path + move, new_state))
    
    return None

# Initial board
board_str = """
+ + + + + + + +
+ - - - - - - +
+ - $ @ - $ - +
+ - X - @ $ - +
+ - - - X * - +
+ + + + + + + +
"""

initial_board = parse_board(board_str)
solution = solve_sokoban(initial_board)
print(solution)