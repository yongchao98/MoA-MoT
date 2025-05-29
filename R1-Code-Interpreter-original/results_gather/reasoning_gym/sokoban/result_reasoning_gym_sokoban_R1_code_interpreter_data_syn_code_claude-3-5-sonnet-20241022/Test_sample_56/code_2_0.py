def parse_board(board_str):
    return [list(row.strip().split()) for row in board_str.strip().split('\n')]

def get_positions(board):
    player = None
    boxes = set()
    goals = set()
    rows, cols = len(board), len(board[0])
    
    for i in range(rows):
        for j in range(cols):
            if board[i][j] == '@':  # Only consider the first '@' as the player
                if player is None:  # Take the first player position
                    player = (i, j)
            elif board[i][j] == '*':  # Player on goal
                if player is None:  # Take the first player position
                    player = (i, j)
                goals.add((i, j))
            elif board[i][j] == '$':  # Box
                boxes.add((i, j))
            elif board[i][j] == 'X':  # Goal
                goals.add((i, j))
    
    return player, boxes, goals

def is_valid_move(board, pos, boxes, direction):
    rows, cols = len(board), len(board[0])
    new_pos = (pos[0] + direction[0], pos[1] + direction[1])
    
    # Check boundaries
    if not (0 <= new_pos[0] < rows and 0 <= new_pos[1] < cols):
        return False, None
    
    # Check if hitting wall
    if board[new_pos[0]][new_pos[1]] == '+':
        return False, None
    
    # If moving into a box
    if new_pos in boxes:
        box_new_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        # Check box new position boundaries
        if not (0 <= box_new_pos[0] < rows and 0 <= box_new_pos[1] < cols):
            return False, None
        # Check if box can be pushed
        if board[box_new_pos[0]][box_new_pos[1]] == '+' or box_new_pos in boxes:
            return False, None
        return True, (new_pos, box_new_pos)
    
    return True, None

def is_deadlock(board, boxes, goals):
    rows, cols = len(board), len(board[0])
    for box in boxes:
        if box in goals:
            continue
        # Check if box is against a wall and not on a goal
        if (board[box[0]][box[1]-1] == '+' and board[box[0]][box[1]+1] == '+') or \
           (board[box[0]-1][box[1]] == '+' and board[box[0]+1][box[1]] == '+'):
            return True
    return False

def solve_sokoban(board):
    player, boxes, goals = get_positions(board)
    if player is None or not boxes or not goals:
        return None
        
    visited = set()
    queue = [(player, frozenset(boxes), "")]
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    while queue:
        curr_pos, curr_boxes, path = queue.pop(0)
        
        if set(curr_boxes) == goals:
            return path
            
        state = (curr_pos, curr_boxes)
        if state in visited:
            continue
            
        visited.add(state)
        
        for move, dir in directions.items():
            valid, box_move = is_valid_move(board, curr_pos, curr_boxes, dir)
            if valid:
                new_pos = (curr_pos[0] + dir[0], curr_pos[1] + dir[1])
                new_boxes = set(curr_boxes)
                
                if box_move:
                    old_box_pos, new_box_pos = box_move
                    new_boxes.remove(old_box_pos)
                    new_boxes.add(new_box_pos)
                    
                    # Skip if deadlock is detected
                    if is_deadlock(board, new_boxes, goals):
                        continue
                
                queue.append((new_pos, frozenset(new_boxes), path + move))
    
    return None

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