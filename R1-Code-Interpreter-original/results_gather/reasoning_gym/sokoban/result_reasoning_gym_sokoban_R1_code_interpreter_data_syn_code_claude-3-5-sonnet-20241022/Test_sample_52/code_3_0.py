def parse_board(board_str):
    return [list(row.strip()) for row in board_str.strip().split('\n')]

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['@', '*']:
                return (i, j)
    return None

def get_boxes_and_goals(board):
    boxes = set()
    goals = set()
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['$', '@']:
                boxes.add((i, j))
            if board[i][j] in ['X', '*']:
                goals.add((i, j))
    return boxes, goals

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def is_deadlock(board, box_pos, goals):
    # Simple deadlock detection: box against wall not on goal
    if box_pos not in goals:
        row, col = box_pos
        if (board[row-1][col] == '+' and board[row][col-1] == '+') or \
           (board[row-1][col] == '+' and board[row][col+1] == '+') or \
           (board[row+1][col] == '+' and board[row][col-1] == '+') or \
           (board[row+1][col] == '+' and board[row][col+1] == '+'):
            return True
    return False

def is_valid_move(board, pos, direction, boxes):
    new_pos = (pos[0] + direction[0], pos[1] + direction[1])
    
    # Check wall collision
    if board[new_pos[0]][new_pos[1]] == '+':
        return False
        
    # Check box pushing
    if new_pos in boxes:
        box_new_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        if (board[box_new_pos[0]][box_new_pos[1]] == '+' or 
            box_new_pos in boxes):
            return False
    return True

def make_move(pos, boxes, direction, goals):
    new_pos = (pos[0] + direction[0], pos[1] + direction[1])
    new_boxes = set(boxes)
    
    if new_pos in boxes:
        box_new_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        new_boxes.remove(new_pos)
        new_boxes.add(box_new_pos)
    
    return new_pos, new_boxes

def solve_sokoban(board):
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    player_pos = get_player_pos(board)
    boxes, goals = get_boxes_and_goals(board)
    
    from collections import deque
    queue = deque([(player_pos, boxes, '')])
    visited = set()
    
    max_moves = 30  # Limit search depth
    
    while queue:
        pos, current_boxes, path = queue.popleft()
        
        if len(path) > max_moves:
            continue
            
        state = (pos, tuple(sorted(current_boxes)))
        if state in visited:
            continue
        visited.add(state)
        
        if current_boxes == goals:
            return path
            
        for move, direction in directions.items():
            if is_valid_move(board, pos, direction, current_boxes):
                new_pos, new_boxes = make_move(pos, current_boxes, direction, goals)
                
                # Skip if any box is in deadlock
                deadlocked = False
                for box in new_boxes:
                    if is_deadlock(board, box, goals):
                        deadlocked = True
                        break
                if deadlocked:
                    continue
                
                queue.append((new_pos, new_boxes, path + move))
    
    return None

# Initial board
board_str = """
+ + + + + + + +
+ - - - - - - +
+ - $ @ - $ - +
+ - X - @ $ - +
+ - - - X * - +
+ + + + + + + +
"""

initial_board = parse_board(board_str)
solution = solve_sokoban(initial_board)
print(solution)