from collections import deque
import copy

def is_valid(board, row, col):
    return 0 <= row < len(board) and 0 <= col < len(board[0]) and board[row][col] != '+'

def get_box_positions(board):
    boxes = set()
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
    return boxes

def get_goals(board):
    goals = set()
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['X', '$', '%']:
                goals.add((i, j))
    return goals

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def is_solved(boxes, goals):
    return boxes == goals

def get_next_states(board, player_pos, boxes):
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    next_states = []
    
    for direction, (dr, dc) in directions.items():
        new_player_row = player_pos[0] + dr
        new_player_col = player_pos[1] + dc
        
        if not is_valid(board, new_player_row, new_player_col):
            continue
            
        new_boxes = boxes.copy()
        if (new_player_row, new_player_col) in boxes:
            new_box_row = new_player_row + dr
            new_box_col = new_player_col + dc
            
            if not is_valid(board, new_box_row, new_box_col) or (new_box_row, new_box_col) in boxes:
                continue
                
            new_boxes.remove((new_player_row, new_player_col))
            new_boxes.add((new_box_row, new_box_col))
        
        next_states.append((direction, (new_player_row, new_player_col), new_boxes))
    
    return next_states

def solve_sokoban(board):
    initial_boxes = get_box_positions(board)
    goals = get_goals(board)
    player_pos = get_player_pos(board)
    
    visited = set()
    queue = deque([(player_pos, initial_boxes, "")])
    
    while queue:
        current_pos, current_boxes, path = queue.popleft()
        
        if is_solved(current_boxes, goals):
            return path
            
        state_key = (current_pos, tuple(sorted(current_boxes)))
        if state_key in visited:
            continue
        visited.add(state_key)
        
        for direction, new_pos, new_boxes in get_next_states(board, current_pos, current_boxes):
            queue.append((new_pos, new_boxes, path + direction))
    
    return None

# Parse the input puzzle
board = []
puzzle = [
    "+ + + + + + + +",
    "+ * @ - - - X +",
    "+ @ $ + - @ - +",
    "+ X + X - - - +",
    "+ + + X - X - +",
    "+ + X - - - - +",
    "+ - - @ X @ - +",
    "+ - @ - - @ - +",
    "+ - - - - - - +",
    "+ + + + + + + +"
]

# Remove spaces from the input
board = [''.join(row.split()) for row in puzzle]

solution = solve_sokoban(board)
print(solution)