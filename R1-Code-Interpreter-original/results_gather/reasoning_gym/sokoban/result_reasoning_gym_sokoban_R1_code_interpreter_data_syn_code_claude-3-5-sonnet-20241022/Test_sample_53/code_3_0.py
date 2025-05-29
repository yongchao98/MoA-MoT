from collections import deque

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def get_boxes(board):
    boxes = set()
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
    return boxes

def get_goals(board):
    goals = set()
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['X', '%', '$']:
                goals.add((i, j))
    return goals

def is_valid_pos(pos, board):
    return (0 <= pos[0] < len(board) and 
            0 <= pos[1] < len(board[0]) and 
            board[pos[0]][pos[1]] != '+')

def make_move(board, player_pos, direction, boxes):
    new_pos = (player_pos[0] + direction[0], player_pos[1] + direction[1])
    
    if not is_valid_pos(new_pos, board):
        return None, None
    
    new_boxes = boxes.copy()
    new_board = [row[:] for row in board]
    
    # Moving to a box
    if new_pos in boxes:
        box_new_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        if not is_valid_pos(box_new_pos, board) or box_new_pos in boxes:
            return None, None
        new_boxes.remove(new_pos)
        new_boxes.add(box_new_pos)
    
    return new_pos, new_boxes

def encode_state(player_pos, boxes):
    boxes_tuple = tuple(sorted(boxes))
    return (player_pos, boxes_tuple)

def solve_sokoban(board):
    player_pos = get_player_pos(board)
    boxes = get_boxes(board)
    goals = get_goals(board)
    
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    queue = deque([(player_pos, boxes, "")])
    seen = {encode_state(player_pos, boxes)}
    
    while queue:
        current_pos, current_boxes, path = queue.popleft()
        
        # Check if solved
        if all(box in goals for box in current_boxes):
            return path
        
        for move, dir in directions.items():
            new_pos, new_boxes = make_move(board, current_pos, dir, current_boxes)
            
            if new_pos is not None:
                state = encode_state(new_pos, new_boxes)
                if state not in seen:
                    seen.add(state)
                    queue.append((new_pos, new_boxes, path + move))
    
    return None

# Initialize the board
board = [
    list("+++++++" ),
    list("+++++X+" ),
    list("++++-@+" ),
    list("+++-@*+" ),
    list("++$-X$+" ),
    list("+-+-+++" ),
    list("+@+++++" ),
    list("+X+++++" ),
    list("+++$+$+" ),
    list("+++++++" )
]

solution = solve_sokoban(board)
print(solution if solution else "No solution found")