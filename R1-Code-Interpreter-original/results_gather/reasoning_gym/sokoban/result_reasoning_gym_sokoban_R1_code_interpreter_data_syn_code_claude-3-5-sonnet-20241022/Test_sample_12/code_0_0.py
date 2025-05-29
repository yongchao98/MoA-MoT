from collections import deque
import heapq
from copy import deepcopy

def parse_board(board_str):
    return [list(row.strip()) for row in board_str.strip().split('\n')]

def find_player(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def find_boxes(board):
    boxes = set()
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
    return boxes

def find_goals(board):
    goals = set()
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['X', '$', '%']:
                goals.add((i, j))
    return goals

def is_valid_move(board, player, direction):
    new_pos = (player[0] + direction[0], player[1] + direction[1])
    if board[new_pos[0]][new_pos[1]] == '+':
        return False
    if board[new_pos[0]][new_pos[1]] in ['@', '$']:
        box_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        if board[box_pos[0]][box_pos[1]] in ['+', '@', '$']:
            return False
    return True

def make_move(board, player, boxes, direction, move_char):
    new_board = [row[:] for row in board]
    new_boxes = set(boxes)
    new_pos = (player[0] + direction[0], player[1] + direction[1])
    
    # Move box if necessary
    if (new_pos[0], new_pos[1]) in boxes:
        box_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        new_boxes.remove(new_pos)
        new_boxes.add(box_pos)
    
    return new_board, new_pos, new_boxes, move_char

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def heuristic(boxes, goals):
    if not boxes or not goals:
        return 0
    total = 0
    unmatched_boxes = list(boxes)
    unmatched_goals = list(goals)
    while unmatched_boxes:
        min_dist = float('inf')
        min_box = None
        min_goal = None
        for box in unmatched_boxes:
            for goal in unmatched_goals:
                dist = manhattan_distance(box, goal)
                if dist < min_dist:
                    min_dist = dist
                    min_box = box
                    min_goal = goal
        total += min_dist
        unmatched_boxes.remove(min_box)
        unmatched_goals.remove(min_goal)
    return total

def solve_sokoban(board_str):
    board = parse_board(board_str)
    player = find_player(board)
    boxes = find_boxes(board)
    goals = find_goals(board)
    
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    start_state = (player, frozenset(boxes))
    visited = set([start_state])
    pq = [(0, '', start_state)]
    
    while pq:
        _, path, (curr_player, curr_boxes) = heapq.heappop(pq)
        
        if curr_boxes == goals:
            return path
        
        for move, dir in directions.items():
            if is_valid_move(board, curr_player, dir):
                _, new_player, new_boxes, _ = make_move(board, curr_player, curr_boxes, dir, move)
                new_state = (new_player, frozenset(new_boxes))
                
                if new_state not in visited:
                    visited.add(new_state)
                    priority = len(path) + 1 + heuristic(new_boxes, goals)
                    heapq.heappush(pq, (priority, path + move, new_state))
    
    return None

# Input puzzle
puzzle = """
+ + + + + + + + + +  
+ $ + + + X - @ - +  
+ + + + + X @ @ - +  
+ X + + X - - X - +  
+ - + X @ - - X - +  
+ @ - - X @ - - - +  
+ - - - @ - @ - - +  
+ - - - - - - * - +  
+ + + + + + + + + +  
"""

solution = solve_sokoban(puzzle)
print(solution)