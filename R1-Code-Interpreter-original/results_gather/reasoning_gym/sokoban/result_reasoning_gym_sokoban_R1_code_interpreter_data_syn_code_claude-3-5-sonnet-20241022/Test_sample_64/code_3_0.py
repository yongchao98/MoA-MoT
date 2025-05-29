from collections import deque

def parse_board(board_str):
    # Split and clean the board string
    board = [list(row.strip().split()) for row in board_str.strip().split('\n')]
    player_pos = None
    boxes = set()
    goals = set()
    
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == '*':  # player on empty
                player_pos = (i, j)
            elif board[i][j] == '%':  # player on goal
                player_pos = (i, j)
                goals.add((i, j))
            elif board[i][j] == '@':  # box
                boxes.add((i, j))
            elif board[i][j] == 'X':  # goal
                goals.add((i, j))
            elif board[i][j] == '$':  # box on goal
                boxes.add((i, j))
                goals.add((i, j))
    
    return board, player_pos, boxes, goals

def is_valid_move(board, player_pos, boxes, direction):
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    dx, dy = moves[direction]
    new_pos = (player_pos[0] + dx, player_pos[1] + dy)
    
    # Check boundaries
    if (new_pos[0] < 0 or new_pos[0] >= len(board) or 
        new_pos[1] < 0 or new_pos[1] >= len(board[0])):
        return False
    
    # Check wall collision
    if board[new_pos[0]][new_pos[1]] == '+':
        return False
    
    # Check box movement
    if new_pos in boxes:
        box_new_pos = (new_pos[0] + dx, new_pos[1] + dy)
        # Check boundaries for box
        if (box_new_pos[0] < 0 or box_new_pos[0] >= len(board) or 
            box_new_pos[1] < 0 or box_new_pos[1] >= len(board[0])):
            return False
        # Check wall or other box collision
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
    max_moves = 200  # Limit the search depth
    
    while queue:
        player_pos, boxes, path = queue.popleft()
        
        if len(path) > max_moves:  # Skip if path is too long
            continue
            
        if boxes == goals:
            return path
            
        state = (player_pos, frozenset(boxes))
        if state in visited:
            continue
            
        visited.add(state)
        
        # Try each direction
        for direction in ['U', 'D', 'L', 'R']:
            if is_valid_move(board, player_pos, boxes, direction):
                new_player, new_boxes = make_move(player_pos, boxes, direction)
                new_path = path + direction
                queue.append((new_player, new_boxes, new_path))
    
    return None

# Test with the puzzle
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