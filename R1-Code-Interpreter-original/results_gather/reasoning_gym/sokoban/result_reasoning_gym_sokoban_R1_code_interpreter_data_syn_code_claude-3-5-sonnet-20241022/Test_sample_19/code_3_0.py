def parse_board(board_str):
    # Remove spaces and create clean board representation
    board = []
    for line in board_str.strip().split('\n'):
        row = ''.join(c for c in line if c != ' ')
        board.append(list(row))
    return board

def get_initial_state(board):
    player = None
    boxes = set()
    goals = set()
    
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['@', '%']:
                player = (i, j)
                if board[i][j] == '%':
                    goals.add((i, j))
            elif board[i][j] == '$':
                boxes.add((i, j))
                goals.add((i, j))
            elif board[i][j] == 'X':
                goals.add((i, j))
    return player, boxes, goals

def is_valid_move(board, pos):
    return (0 <= pos[0] < len(board) and 
            0 <= pos[1] < len(board[0]) and 
            board[pos[0]][pos[1]] != '+')

def get_next_states(pos, boxes, board):
    moves = []
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    for move, (dy, dx) in directions.items():
        new_pos = (pos[0] + dy, pos[1] + dx)
        
        if not is_valid_move(board, new_pos):
            continue
            
        if new_pos in boxes:
            # Try to push the box
            box_new_pos = (new_pos[0] + dy, new_pos[1] + dx)
            if not is_valid_move(board, box_new_pos) or box_new_pos in boxes:
                continue
                
            # Create new box configuration
            new_boxes = set(boxes)
            new_boxes.remove(new_pos)
            new_boxes.add(box_new_pos)
            moves.append((new_pos, new_boxes, move))
        else:
            # Simple player move
            moves.append((new_pos, boxes, move))
    
    return moves

def solve_sokoban(board_str):
    from collections import deque
    
    # Initialize the board and state
    board = parse_board(board_str)
    player, boxes, goals = get_initial_state(board)
    
    # BFS queue and visited states
    queue = deque([(player, boxes, "")])
    visited = set()
    
    while queue:
        current_pos, current_boxes, path = queue.popleft()
        
        # Check if we reached the goal state
        if current_boxes == goals:
            return path
        
        # Generate state key for visited check
        state = (current_pos, frozenset(current_boxes))
        if state in visited:
            continue
        visited.add(state)
        
        # Try all possible moves
        for next_pos, next_boxes, move in get_next_states(current_pos, current_boxes, board):
            next_state = (next_pos, frozenset(next_boxes))
            if next_state not in visited:
                queue.append((next_pos, next_boxes, path + move))
    
    return None

# Test the solution
puzzle = """+ + + + + + +
+ $ + - - X +
+ - - - - - +
+ - - - @ - +
+ - X - @ - +
+ - - - X - +
+ - - @ X X +
+ - - - @ @ +
+ - - - @ % +
+ + + + + + +"""

solution = solve_sokoban(puzzle)
print(solution)