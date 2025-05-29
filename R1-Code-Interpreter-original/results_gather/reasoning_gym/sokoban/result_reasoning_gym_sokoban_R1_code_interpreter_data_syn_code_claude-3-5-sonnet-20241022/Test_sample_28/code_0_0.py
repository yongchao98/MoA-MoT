from collections import deque
import copy

def is_valid_pos(x, y, width, height):
    return 0 <= x < height and 0 <= y < width

def is_deadlock(board, box_x, box_y):
    # Simple corner deadlock detection
    if (board[box_x-1][box_y] == '+' and board[box_x][box_y-1] == '+') or \
       (board[box_x-1][box_y] == '+' and board[box_x][box_y+1] == '+') or \
       (board[box_x+1][box_y] == '+' and board[box_x][box_y-1] == '+') or \
       (board[box_x+1][box_y] == '+' and board[box_x][box_y+1] == '+'):
        return True
    return False

def get_state_key(players, boxes):
    return (tuple(sorted(players)), tuple(sorted(boxes)))

def solve_sokoban(board):
    height = len(board)
    width = len(board[0])
    
    # Find initial positions
    players = set()
    boxes = set()
    goals = set()
    
    for i in range(height):
        for j in range(width):
            if board[i][j] in ['*', '%']:
                players.add((i, j))
            elif board[i][j] in ['@', '$']:
                boxes.add((i, j))
            if board[i][j] in ['X', '$', '%']:
                goals.add((i, j))
    
    # BFS queue: (players, boxes, moves)
    queue = deque([(players, boxes, "")])
    visited = {get_state_key(players, boxes)}
    
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    while queue:
        current_players, current_boxes, moves = queue.popleft()
        
        # Check if won
        if all(box in goals for box in current_boxes):
            return moves
        
        # Try each player's moves
        for player in current_players:
            px, py = player
            
            for move, (dx, dy) in directions.items():
                new_px, new_py = px + dx, py + dy
                
                if not is_valid_pos(new_px, new_py, width, height) or board[new_px][new_py] == '+':
                    continue
                
                new_players = set(current_players)
                new_players.remove(player)
                new_players.add((new_px, new_py))
                new_boxes = set(current_boxes)
                
                # If pushing a box
                if (new_px, new_py) in current_boxes:
                    box_new_x, box_new_y = new_px + dx, new_py + dy
                    
                    if not is_valid_pos(box_new_x, box_new_y, width, height) or \
                       board[box_new_x][box_new_y] == '+' or \
                       (box_new_x, box_new_y) in current_boxes or \
                       is_deadlock(board, box_new_x, box_new_y):
                        continue
                    
                    new_boxes.remove((new_px, new_py))
                    new_boxes.add((box_new_x, box_new_y))
                
                new_state = get_state_key(new_players, new_boxes)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_players, new_boxes, moves + move))
    
    return None

# Initialize the board
board = [
    list("++++++++"),
    list("+ --*X +"),
    list("++$-@- +"),
    list("+$X--@-+"),
    list("+X--X- +"),
    list("+-@--@-+"),
    list("+------+"),
    list("++++++++")
]

solution = solve_sokoban(board)
print(solution)