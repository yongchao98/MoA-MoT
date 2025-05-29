from collections import deque
import copy

def parse_board(board_str):
    board = [list(row.strip()) for row in board_str.strip().split('\n')]
    player_pos = None
    boxes = set()
    goals = set()
    
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['@', '*']:
                player_pos = (i, j)
                if board[i][j] == '*':
                    goals.add((i, j))
            elif board[i][j] == '@':
                boxes.add((i, j))
            elif board[i][j] == '$':
                boxes.add((i, j))
                goals.add((i, j))
            elif board[i][j] == 'X':
                goals.add((i, j))
    
    return board, player_pos, boxes, goals

def is_valid_move(board, player_pos, boxes, direction):
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    dx, dy = moves[direction]
    new_pos = (player_pos[0] + dx, player_pos[1] + dy)
    
    # Check if player can move
    if board[new_pos[0]][new_pos[1]] == '+':
        return False
        
    # If moving into box, check if box can move
    if new_pos in boxes:
        box_new_pos = (new_pos[0] + dx, new_pos[1] + dy)
        if board[box_new_pos[0]][box_new_pos[1]] == '+' or box_new_pos in boxes:
            return False
    return True

def make_move(player_pos, boxes, direction):
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    dx, dy = moves[direction]
    new_pos = (player_pos[0] + dx, player_pos[1] + dy)
    new_boxes = set(boxes)
    
    if new_pos in boxes:
        box_new_pos = (new_pos[0] + dx, new_pos[1] + dy)
        new_boxes.remove(new_pos)
        new_boxes.add(box_new_pos)
    
    return new_pos, new_boxes

def solve_sokoban(board_str):
    board, initial_player, initial_boxes, goals = parse_board(board_str)
    visited = set()
    queue = deque([(initial_player, initial_boxes, "")])
    
    while queue:
        player_pos, boxes, path = queue.popleft()
        state = (player_pos, frozenset(boxes))
        
        if boxes == goals:
            return path
            
        if state in visited:
            continue
            
        visited.add(state)
        
        for direction in ['U', 'D', 'L', 'R']:
            if is_valid_move(board, player_pos, boxes, direction):
                new_player, new_boxes = make_move(player_pos, boxes, direction)
                queue.append((new_player, new_boxes, path + direction))
    
    return None

# Test the puzzle
puzzle = """
+ + + + + + + +
+ + - X - - - +
+ X - @ X @ * +
+ X @ - - - - +
+ + - - - @ - +
+ - - - - @ - +
+ X @ X - - - +
+ + + + + + + +
"""

solution = solve_sokoban(puzzle)
print(f"<<<{solution}>>>")