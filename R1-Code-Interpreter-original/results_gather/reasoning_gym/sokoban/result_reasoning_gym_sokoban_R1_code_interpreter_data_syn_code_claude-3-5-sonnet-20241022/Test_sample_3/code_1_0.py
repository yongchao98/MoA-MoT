from collections import deque
import copy

def is_valid_pos(x, y, width, height):
    return 0 <= x < height and 0 <= y < width

def get_next_pos(x, y, direction):
    if direction == 'U': return (x-1, y)
    if direction == 'D': return (x+1, y)
    if direction == 'L': return (x, y-1)
    if direction == 'R': return (x, y+1)

def is_box(char):
    return char in ['@', '$']

def is_player(char):
    return char in ['@', '%', '*']

def is_goal(char):
    return char in ['X', '$', '%']

def solve_sokoban(board):
    height = len(board)
    width = len(board[0])
    
    # Find initial positions
    player_pos = None
    boxes = set()
    goals = set()
    
    for i in range(height):
        for j in range(width):
            if is_player(board[i][j]):
                player_pos = (i, j)
            if is_box(board[i][j]):
                boxes.add((i, j))
            if is_goal(board[i][j]):
                goals.add((i, j))
    
    # BFS
    queue = deque([(player_pos, boxes, "")])
    visited = set()
    
    while queue:
        curr_player, curr_boxes, path = queue.popleft()
        
        # Check if solved
        if all(box in goals for box in curr_boxes):
            return path
        
        # Try all directions
        for direction in ['U', 'D', 'L', 'R']:
            next_player = get_next_pos(curr_player[0], curr_player[1], direction)
            
            if not is_valid_pos(next_player[0], next_player[1], width, height):
                continue
                
            if board[next_player[0]][next_player[1]] == '+':
                continue
            
            new_boxes = set(curr_boxes)
            
            # If pushing a box
            if next_player in curr_boxes:
                next_box = get_next_pos(next_player[0], next_player[1], direction)
                
                if not is_valid_pos(next_box[0], next_box[1], width, height):
                    continue
                    
                if board[next_box[0]][next_box[1]] == '+':
                    continue
                    
                if next_box in curr_boxes:
                    continue
                    
                new_boxes.remove(next_player)
                new_boxes.add(next_box)
            
            state = (next_player, tuple(sorted(new_boxes)))
            if state not in visited:
                visited.add(state)
                queue.append((next_player, new_boxes, path + direction))
    
    return None

# Initial board (fixed representation)
board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', '-', '-', '-', '-', '-', '+', '$', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '+', '+', '+'],
    ['+', 'X', '-', 'X', '-', '-', '-', '$', '+', '+'],
    ['+', '-', '@', '-', '-', '-', '-', '+', '+', '+'],
    ['+', 'X', '-', '@', '@', '@', '+', '+', '+', '+'],
    ['+', '+', 'X', '-', '@', '%', '$', '+', '+', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(board)
print(solution)