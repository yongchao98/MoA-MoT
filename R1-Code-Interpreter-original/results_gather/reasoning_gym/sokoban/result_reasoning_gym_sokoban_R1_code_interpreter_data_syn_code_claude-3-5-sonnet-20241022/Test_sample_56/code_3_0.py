def parse_board(board_str):
    lines = [line.strip().split() for line in board_str.strip().split('\n')]
    return lines

def get_initial_state(board):
    rows, cols = len(board), len(board[0])
    player = None
    boxes = set()
    goals = set()
    
    # Find player, boxes, and goals
    for i in range(rows):
        for j in range(cols):
            cell = board[i][j]
            if cell == '@':  # Player
                player = (i, j)
            elif cell == '*':  # Player on goal
                player = (i, j)
                goals.add((i, j))
            elif cell == '$':  # Box
                boxes.add((i, j))
            elif cell == 'X':  # Goal
                goals.add((i, j))
    
    return player, boxes, goals

def is_valid_push(board, player_pos, box_pos, direction):
    new_box_pos = (box_pos[0] + direction[0], box_pos[1] + direction[1])
    rows, cols = len(board), len(board[0])
    
    # Check boundaries
    if not (0 <= new_box_pos[0] < rows and 0 <= new_box_pos[1] < cols):
        return False
    
    # Check if new box position is valid
    if board[new_box_pos[0]][new_box_pos[1]] == '+':
        return False
        
    return True

def get_next_states(board, player_pos, boxes, path):
    next_states = []
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    for move, (dy, dx) in directions.items():
        new_pos = (player_pos[0] + dy, player_pos[1] + dx)
        
        # Check if new position is within bounds
        if not (0 <= new_pos[0] < len(board) and 0 <= new_pos[1] < len(board[0])):
            continue
            
        # Check if new position hits a wall
        if board[new_pos[0]][new_pos[1]] == '+':
            continue
            
        new_boxes = set(boxes)
        if new_pos in boxes:  # If moving into a box
            if is_valid_push(board, player_pos, new_pos, (dy, dx)):
                new_boxes.remove(new_pos)
                new_boxes.add((new_pos[0] + dy, new_pos[1] + dx))
                next_states.append((new_pos, new_boxes, path + move))
        else:  # If moving to empty space
            next_states.append((new_pos, new_boxes, path + move))
            
    return next_states

def solve_sokoban(board):
    from collections import deque
    
    player, boxes, goals = get_initial_state(board)
    if not player or not boxes or not goals:
        return None
        
    queue = deque([(player, boxes, "")])
    visited = set()
    
    while queue:
        curr_pos, curr_boxes, path = queue.popleft()
        
        # Check if we reached the goal
        if curr_boxes == goals:
            return path
            
        # Generate state key
        state = (curr_pos, frozenset(curr_boxes))
        if state in visited:
            continue
            
        visited.add(state)
        
        # Get next possible states
        next_states = get_next_states(board, curr_pos, curr_boxes, path)
        for next_state in next_states:
            queue.append(next_state)
    
    return None

# Test the puzzle
puzzle = """
+ + + + + + + + +
+ - - - + X - $ +
+ X - - + $ - * +
+ - - @ - - X - +
+ @ $ - - - @ - +
+ - - - - - - - +
+ + + + + + + + +
"""

solution = solve_sokoban(parse_board(puzzle))
print(solution)