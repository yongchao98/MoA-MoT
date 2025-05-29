from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_positions(board, pos, boxes):
    rows, cols = len(board), len(board[0])
    moves = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    next_positions = []
    
    for dx, dy, direction in moves:
        new_x, new_y = pos[0] + dx, pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or board[new_x][new_y] == '+':
            continue
            
        if (new_x, new_y) in boxes:
            box_new_x, box_new_y = new_x + dx, new_y + dy
            if (is_valid(box_new_x, box_new_y, rows, cols) and 
                board[box_new_x][box_new_y] != '+' and 
                (box_new_x, box_new_y) not in boxes):
                next_positions.append((new_x, new_y, direction, (new_x, new_y), (box_new_x, box_new_y)))
        else:
            next_positions.append((new_x, new_y, direction, None, None))
            
    return next_positions

def solve_sokoban(board):
    rows, cols = len(board), len(board[0])
    players = []
    boxes = set()
    goals = set()
    
    # Parse the board
    for i in range(rows):
        for j in range(cols):
            if board[i][j] in ['*', '%']:
                players.append((i, j))
            elif board[i][j] == '@':
                boxes.add((i, j))
            elif board[i][j] == '$':
                boxes.add((i, j))
                goals.add((i, j))
            elif board[i][j] == 'X':
                goals.add((i, j))
    
    # BFS
    queue = deque([(players, boxes, "")])
    visited = set()
    
    while queue:
        current_players, current_boxes, path = queue.popleft()
        
        # Check if all boxes are on goals
        if all(box in goals for box in current_boxes):
            return path
        
        state = (tuple(sorted(current_players)), tuple(sorted(current_boxes)))
        if state in visited:
            continue
        visited.add(state)
        
        # Try moving each player
        for player_idx, player_pos in enumerate(current_players):
            next_positions = get_next_positions(board, player_pos, current_boxes)
            
            for new_x, new_y, direction, box_pos, new_box_pos in next_positions:
                new_players = list(current_players)
                new_players[player_idx] = (new_x, new_y)
                new_boxes = set(current_boxes)
                
                if box_pos:
                    new_boxes.remove(box_pos)
                    new_boxes.add(new_box_pos)
                
                queue.append((new_players, new_boxes, path + direction))
    
    return "No solution found"

# Initialize the board
board = [
    ['+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '+'],
    ['+', '@', '-', '$', 'X', '+'],
    ['+', '-', '-', '@', '-', '+'],
    ['+', '-', 'X', '*', '-', '+'],
    ['+', '-', '-', '-', '-', '+'],
    ['+', 'X', '-', '@', '-', '+'],
    ['+', '$', '-', '-', '-', '+'],
    ['+', 'X', '@', '-', '+', '+'],
    ['+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(board)
print(solution)