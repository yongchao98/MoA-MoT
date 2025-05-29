from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_pos(x, y, direction):
    if direction == 'U': return (x-1, y)
    if direction == 'D': return (x+1, y)
    if direction == 'L': return (x, y-1)
    if direction == 'R': return (x, y+1)

def convert_to_state(board):
    rows, cols = len(board), len(board[0])
    player_pos = None
    boxes = set()
    goals = set()
    
    for i in range(rows):
        for j in range(cols):
            if board[i][j] in ['*', '%', '@']:
                if player_pos is None:  # First player position found
                    player_pos = (i, j)
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
            if board[i][j] in ['X', '%', '$']:
                goals.add((i, j))
    
    return player_pos, boxes, goals

def solve_sokoban(board):
    initial_board = [list(row) for row in board.split('\n')]
    rows, cols = len(initial_board), len(initial_board[0])
    
    player_pos, boxes, goals = convert_to_state(initial_board)
    
    # BFS queue: (player_pos, boxes, path)
    queue = deque([(player_pos, frozenset(boxes), "")])
    visited = set()
    
    while queue:
        curr_player, curr_boxes, path = queue.popleft()
        
        if all(box in goals for box in curr_boxes):
            return path
        
        state = (curr_player, frozenset(curr_boxes))
        if state in visited:
            continue
        visited.add(state)
        
        for direction in ['U', 'D', 'L', 'R']:
            next_player = get_next_pos(*curr_player, direction)
            
            if not is_valid(*next_player, rows, cols):
                continue
                
            if initial_board[next_player[0]][next_player[1]] == '+':
                continue
                
            new_boxes = set(curr_boxes)
            
            # If there's a box in the next position
            if next_player in curr_boxes:
                next_box = get_next_pos(*next_player, direction)
                
                if not is_valid(*next_box, rows, cols):
                    continue
                    
                if initial_board[next_box[0]][next_box[1]] == '+':
                    continue
                    
                if next_box in curr_boxes:
                    continue
                    
                new_boxes.remove(next_player)
                new_boxes.add(next_box)
            
            queue.append((next_player, frozenset(new_boxes), path + direction))
    
    return None

# Input puzzle
puzzle = """+ + + + + + + +
+ + + + $ + + +
+ + X + + $ + +
+ X - $ + + + +
+ - @ + + + + +
+ - - - $ + + +
+ @ - - - - + +
+ * @ - X + + +
+ + + + + + + +"""

solution = solve_sokoban(puzzle)
print(solution)