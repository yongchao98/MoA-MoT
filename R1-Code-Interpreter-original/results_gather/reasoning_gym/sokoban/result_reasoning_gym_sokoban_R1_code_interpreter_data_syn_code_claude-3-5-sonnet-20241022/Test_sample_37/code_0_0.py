from collections import deque
import copy

def is_valid(board, r, c):
    return 0 <= r < len(board) and 0 <= c < len(board[0]) and board[r][c] != '+'

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def get_boxes(board):
    boxes = set()
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
    return boxes

def get_goals(board):
    goals = set()
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['X', '$', '%']:
                goals.add((i, j))
    return goals

def is_win(boxes, goals):
    return boxes == goals

def solve_sokoban(initial_board):
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    player_pos = get_player_pos(initial_board)
    boxes = get_boxes(initial_board)
    goals = get_goals(initial_board)
    
    queue = deque([(player_pos, boxes, "")])
    visited = set()
    
    while queue:
        pos, current_boxes, path = queue.popleft()
        
        if is_win(current_boxes, goals):
            return path
            
        state = (pos, tuple(sorted(current_boxes)))
        if state in visited:
            continue
        visited.add(state)
        
        for move, (dr, dc) in directions.items():
            new_r, new_c = pos[0] + dr, pos[1] + dc
            
            if not is_valid(initial_board, new_r, new_c):
                continue
                
            if (new_r, new_c) in current_boxes:
                box_new_r, box_new_c = new_r + dr, new_c + dc
                if not is_valid(initial_board, box_new_r, box_new_c):
                    continue
                if (box_new_r, box_new_c) in current_boxes:
                    continue
                    
                new_boxes = set(current_boxes)
                new_boxes.remove((new_r, new_c))
                new_boxes.add((box_new_r, box_new_c))
                queue.append(((new_r, new_c), new_boxes, path + move))
            else:
                queue.append(((new_r, new_c), current_boxes, path + move))
    
    return None

# Initialize the board
board = [
    ['+', '+', '+', '+', '+', '+'],
    ['+', '*', '-', '-', '-', '+'],
    ['+', '-', '@', '-', '$', '+'],
    ['+', '-', '@', 'X', 'X', '+'],
    ['+', '-', '@', '-', 'X', '+'],
    ['+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(board)
print(solution)