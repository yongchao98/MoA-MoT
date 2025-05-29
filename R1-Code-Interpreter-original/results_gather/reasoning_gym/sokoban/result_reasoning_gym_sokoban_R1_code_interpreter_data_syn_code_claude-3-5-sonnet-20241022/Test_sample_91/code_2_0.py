from collections import deque

def is_valid(x, y, board):
    return 0 <= x < len(board) and 0 <= y < len(board[0]) and board[x][y] != '+'

def get_state_string(board, boxes):
    state = []
    for row in board:
        state.append(''.join(row))
    return '\n'.join(state) + '\n' + ','.join(f"{x},{y}" for x, y in sorted(boxes))

def solve_sokoban(initial):
    board = [list(row) for row in initial]
    rows, cols = len(board), len(board[0])
    
    # Find initial positions
    players = set()
    boxes = set()
    goals = set()
    
    for i in range(rows):
        for j in range(cols):
            if board[i][j] in ['@', '*']:
                players.add((i, j))
                if board[i][j] == '*':
                    goals.add((i, j))
            elif board[i][j] == '$':
                boxes.add((i, j))
            elif board[i][j] == 'X':
                goals.add((i, j))
    
    # Directions: Up, Down, Left, Right
    directions = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    
    visited = set()
    queue = deque([(players, boxes, "")])
    
    while queue:
        current_players, current_boxes, moves = queue.popleft()
        
        # Check if won
        if all(box in goals for box in current_boxes):
            return moves
        
        state = get_state_string(board, current_boxes)
        if state in visited:
            continue
        visited.add(state)
        
        # Try each player's moves
        for player in current_players:
            px, py = player
            
            for dx, dy, move in directions:
                new_x, new_y = px + dx, py + dy
                
                if not is_valid(new_x, new_y, board):
                    continue
                
                new_players = set(current_players)
                new_boxes = set(current_boxes)
                
                # If pushing a box
                if (new_x, new_y) in current_boxes:
                    box_x, box_y = new_x + dx, new_y + dy
                    
                    # Check if box can be pushed
                    if not is_valid(box_x, box_y, board) or (box_x, box_y) in current_boxes:
                        continue
                    
                    new_boxes.remove((new_x, new_y))
                    new_boxes.add((box_x, box_y))
                
                # Update player position
                new_players.remove((px, py))
                new_players.add((new_x, new_y))
                
                queue.append((new_players, new_boxes, moves + move))
    
    return None

initial = [
    "++++++",
    "+-@-X+",
    "+-$-++",
    "+X@X-+",
    "+-@--+",
    "+-@X*+",
    "+----+",
    "++++++"
]

solution = solve_sokoban(initial)
print(solution)