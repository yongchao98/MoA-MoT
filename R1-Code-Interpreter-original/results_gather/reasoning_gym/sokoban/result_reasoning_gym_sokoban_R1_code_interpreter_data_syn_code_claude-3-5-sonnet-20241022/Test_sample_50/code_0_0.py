from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_player_positions(board):
    positions = []
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%', '@']:
                positions.append((i, j))
    return positions

def get_box_positions(board):
    positions = []
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['$', '@']:
                positions.append((i, j))
    return positions

def get_goal_positions(board):
    positions = []
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['X', '$']:
                positions.append((i, j))
    return positions

def is_win(boxes, goals):
    return all(box in goals for box in boxes)

def get_next_states(board, player, boxes):
    moves = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    
    for move, dx, dy in moves:
        new_player = (player[0] + dx, player[1] + dy)
        
        if not is_valid(new_player[0], new_player[1], len(board), len(board[0])):
            continue
            
        if board[new_player[0]][new_player[1]] == '+':
            continue
            
        if new_player in boxes:
            new_box = (new_player[0] + dx, new_player[1] + dy)
            if not is_valid(new_box[0], new_box[1], len(board), len(board[0])):
                continue
            if board[new_box[0]][new_box[1]] == '+':
                continue
            if new_box in boxes:
                continue
                
            new_boxes = tuple(new_box if b == new_player else b for b in boxes)
            next_states.append((move, new_player, new_boxes))
        else:
            next_states.append((move, new_player, boxes))
            
    return next_states

def solve_sokoban(board):
    initial_board = [list(row) for row in board.strip().split('\n')]
    players = get_player_positions(initial_board)
    boxes = tuple(get_box_positions(initial_board))
    goals = tuple(get_goal_positions(initial_board))
    
    queue = deque([(players[0], boxes, '')])
    visited = set()
    
    while queue:
        player, current_boxes, path = queue.popleft()
        
        if is_win(current_boxes, goals):
            return path
            
        state = (player, current_boxes)
        if state in visited:
            continue
        visited.add(state)
        
        for move, new_player, new_boxes in get_next_states(initial_board, player, current_boxes):
            queue.append((new_player, new_boxes, path + move))
    
    return None

# Define the puzzle
puzzle = """+ + + + + + +
+ - - - - - +
+ $ @ X - - +
+ - - - @ - +
+ - - - $ - +
+ X - - - - +
+ + @ X @ - +
+ + X - - * +
+ + + + + + +"""

solution = solve_sokoban(puzzle)
print(solution)