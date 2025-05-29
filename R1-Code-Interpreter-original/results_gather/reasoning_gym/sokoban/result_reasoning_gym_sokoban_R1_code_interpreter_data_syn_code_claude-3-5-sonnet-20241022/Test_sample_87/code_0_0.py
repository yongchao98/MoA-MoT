from collections import deque
import copy

def is_valid_pos(x, y, width, height):
    return 0 <= x < height and 0 <= y < width

def get_next_pos(x, y, direction):
    if direction == 'U': return (x-1, y)
    if direction == 'D': return (x+1, y)
    if direction == 'L': return (x, y-1)
    if direction == 'R': return (x, y+1)

def convert_to_state(board):
    player = None
    boxes = set()
    goals = set()
    height = len(board)
    width = len(board[0])
    
    for i in range(height):
        for j in range(width):
            if board[i][j] == '*':
                player = (i, j)
            elif board[i][j] == '@':
                boxes.add((i, j))
            elif board[i][j] == 'X':
                goals.add((i, j))
            elif board[i][j] == '$':
                boxes.add((i, j))
                goals.add((i, j))
            elif board[i][j] == '%':
                player = (i, j)
                goals.add((i, j))
    
    return player, boxes, goals

def is_valid_move(board, player, new_player, boxes, direction):
    if not is_valid_pos(new_player[0], new_player[1], len(board), len(board[0])):
        return False
    if board[new_player[0]][new_player[1]] == '+':
        return False
    
    if new_player in boxes:
        box_new_pos = get_next_pos(new_player[0], new_player[1], direction)
        if not is_valid_pos(box_new_pos[0], box_new_pos[1], len(board), len(board[0])):
            return False
        if board[box_new_pos[0]][box_new_pos[1]] == '+':
            return False
        if box_new_pos in boxes:
            return False
    return True

def solve_sokoban(board):
    player, boxes, goals = convert_to_state(board)
    queue = deque([(player, boxes, "")])
    visited = set()
    
    while queue:
        current_player, current_boxes, path = queue.popleft()
        
        if current_boxes == goals:
            return path
        
        state = (current_player, frozenset(current_boxes))
        if state in visited:
            continue
        visited.add(state)
        
        for direction in ['U', 'D', 'L', 'R']:
            new_player = get_next_pos(current_player[0], current_player[1], direction)
            
            if is_valid_move(board, current_player, new_player, current_boxes, direction):
                new_boxes = set(current_boxes)
                
                if new_player in new_boxes:
                    box_new_pos = get_next_pos(new_player[0], new_player[1], direction)
                    new_boxes.remove(new_player)
                    new_boxes.add(box_new_pos)
                
                queue.append((new_player, new_boxes, path + direction))
    
    return None

# Initialize the board
board = [
    ['+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', 'X', '@', '-', '-', '-', '-', '+'],
    ['+', '+', 'X', '-', '@', '@', '@', '+'],
    ['+', '+', '-', 'X', '*', '-', '-', '+'],
    ['+', '-', '-', '@', '-', '-', 'X', '+'],
    ['+', '-', '-', '-', 'X', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(board)
print(solution)