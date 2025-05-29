def encode_state(board):
    # Encode board state more efficiently
    boxes = set()
    goals = set()
    player = None
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
            if board[i][j] in ['X', '$', '%']:
                goals.add((i, j))
            if board[i][j] in ['*', '%']:
                player = (i, j)
    return (player, frozenset(boxes), frozenset(goals))

def is_deadlock(board, box_pos):
    r, c = box_pos
    # Check if box is in corner
    if (board[r-1][c] == '+' and board[r][c-1] == '+') or \
       (board[r-1][c] == '+' and board[r][c+1] == '+') or \
       (board[r+1][c] == '+' and board[r][c-1] == '+') or \
       (board[r+1][c] == '+' and board[r][c+1] == '+'):
        return True
    return False

def solve_sokoban_optimized(board):
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    moves = {'U': 'U', 'D': 'D', 'L': 'L', 'R': 'R'}
    
    # Find initial state
    initial_state = encode_state(board)
    player_pos = initial_state[0]
    boxes = initial_state[1]
    goals = initial_state[2]
    
    from collections import deque
    queue = deque([(player_pos, boxes, "")])
    visited = set()
    
    while queue:
        pos, current_boxes, path = queue.popleft()
        
        # Check if solved
        if current_boxes == goals:
            return path
            
        # Generate state key
        state = (pos, frozenset(current_boxes))
        if state in visited:
            continue
        visited.add(state)
        
        # Try each direction
        for move, (dr, dc) in directions.items():
            new_r, new_c = pos[0] + dr, pos[1] + dc
            
            # Check if move is valid
            if board[new_r][new_c] == '+':
                continue
                
            new_boxes = set(current_boxes)
            if (new_r, new_c) in current_boxes:
                # Pushing a box
                next_r, next_c = new_r + dr, new_c + dc
                if board[next_r][next_c] != '+' and (next_r, next_c) not in current_boxes:
                    if not is_deadlock(board, (next_r, next_c)):
                        new_boxes.remove((new_r, new_c))
                        new_boxes.add((next_r, next_c))
                        queue.append(((new_r, new_c), frozenset(new_boxes), path + move))
            else:
                # Just moving player
                queue.append(((new_r, new_c), frozenset(new_boxes), path + move))
    
    return None

# Initialize the board
board = [
    ['+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '@', '-', '-', '+'],
    ['+', '-', '-', '-', '-', 'X', '+'],
    ['+', 'X', '@', '-', 'X', '-', '+'],
    ['+', '-', '@', '*', '@', '-', '+'],
    ['+', '-', '-', 'X', '@', 'X', '+'],
    ['+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban_optimized(board)
print(solution)