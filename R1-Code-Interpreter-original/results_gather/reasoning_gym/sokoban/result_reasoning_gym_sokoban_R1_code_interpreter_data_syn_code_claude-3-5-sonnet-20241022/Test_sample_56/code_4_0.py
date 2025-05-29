def get_initial_state(board):
    rows = len(board)
    cols = len(board[0])
    player = None
    boxes = set()
    goals = set()
    
    # Find first player position (ignore others), boxes, and goals
    for i in range(rows):
        for j in range(cols):
            if board[i][j] == '@' and player is None:
                player = (i, j)
            elif board[i][j] == '$':
                boxes.add((i, j))
            elif board[i][j] == 'X':
                goals.add((i, j))
            elif board[i][j] == '*':  # Player on goal
                if player is None:
                    player = (i, j)
                goals.add((i, j))
    
    return player, boxes, goals

def is_valid_move(pos, board):
    return (0 <= pos[0] < len(board) and 
            0 <= pos[1] < len(board[0]) and 
            board[pos[0]][pos[1]] != '+')

def solve_sokoban(board_str):
    # Parse board
    board = [list(row.split()) for row in board_str.strip().split('\n')]
    player, boxes, goals = get_initial_state(board)
    
    # BFS queue and visited states
    from collections import deque
    queue = deque([(player, frozenset(boxes), "")])
    visited = set()
    
    # Directions
    moves = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    while queue:
        curr_pos, curr_boxes, path = queue.popleft()
        
        # Check if solved
        if curr_boxes == goals:
            return path
        
        # Generate state key
        state = (curr_pos, curr_boxes)
        if state in visited:
            continue
        visited.add(state)
        
        # Try each direction
        for move, (dy, dx) in moves.items():
            new_pos = (curr_pos[0] + dy, curr_pos[1] + dx)
            
            # Check if new position is valid
            if not is_valid_move(new_pos, board):
                continue
                
            new_boxes = set(curr_boxes)
            
            # If moving into a box
            if new_pos in curr_boxes:
                box_new_pos = (new_pos[0] + dy, new_pos[1] + dx)
                
                # Check if box can be pushed
                if not is_valid_move(box_new_pos, board) or box_new_pos in curr_boxes:
                    continue
                    
                new_boxes.remove(new_pos)
                new_boxes.add(box_new_pos)
            
            queue.append((new_pos, frozenset(new_boxes), path + move))
    
    return None

# Test puzzle
puzzle = """
+ + + + + + + + +
+ - - - + X - $ +
+ X - - + $ - * +
+ - - @ - - X - +
+ @ $ - - - @ - +
+ - - - - - - - +
+ + + + + + + + +
"""

solution = solve_sokoban(puzzle)
print(solution)