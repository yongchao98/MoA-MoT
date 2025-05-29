from collections import deque

def get_state(board):
    player = None
    boxes = set()
    goals = set()
    walls = set()
    
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['*', '%']:
                player = (i, j)
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
            if board[i][j] in ['X', '$', '%']:
                goals.add((i, j))
            if board[i][j] == '+':
                walls.add((i, j))
    
    return player, boxes, goals, walls

def is_deadlock(pos, boxes, goals, walls, board_size):
    if pos in goals:
        return False
        
    rows, cols = board_size
    x, y = pos
    
    # Corner deadlock
    corners = [
        ((x-1, y), (x, y-1)),
        ((x-1, y), (x, y+1)),
        ((x+1, y), (x, y-1)),
        ((x+1, y), (x, y+1))
    ]
    
    for (wx, wy), (wz, wk) in corners:
        if (wx, wy) in walls and (wz, wk) in walls:
            return True
    
    return False

def make_move(player, boxes, dx, dy, walls, board_size):
    px, py = player
    new_x, new_y = px + dx, py + dy
    
    if (new_x, new_y) in walls:
        return None, None
        
    new_boxes = set(boxes)
    
    if (new_x, new_y) in boxes:
        push_x, push_y = new_x + dx, new_y + dy
        if (push_x, push_y) in walls or (push_x, push_y) in boxes:
            return None, None
        new_boxes.remove((new_x, new_y))
        new_boxes.add((push_x, push_y))
        
    return (new_x, new_y), new_boxes

def solve_sokoban(board):
    player, boxes, goals, walls = get_state(board)
    board_size = (len(board), len(board[0]))
    
    initial_state = (player, frozenset(boxes))
    visited = {initial_state}
    queue = deque([(initial_state, "")])
    
    directions = [
        (0, 1, 'R'),
        (0, -1, 'L'),
        (-1, 0, 'U'),
        (1, 0, 'D')
    ]
    
    max_moves = 100
    
    while queue:
        state, path = queue.popleft()
        
        if len(path) >= max_moves:
            continue
            
        current_player, current_boxes = state
        
        if all(box in goals for box in current_boxes):
            return path
            
        for dx, dy, move in directions:
            new_player, new_boxes = make_move(current_player, current_boxes, dx, dy, walls, board_size)
            
            if new_player is None:
                continue
                
            # Check for deadlocks
            deadlocked = False
            for box in new_boxes:
                if is_deadlock(box, new_boxes, goals, walls, board_size):
                    deadlocked = True
                    break
                    
            if deadlocked:
                continue
                
            new_state = (new_player, frozenset(new_boxes))
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + move))
    
    return None

# Initialize the board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '+', 'X', '+', '$', 'X', '+'],
    ['+', '-', 'X', 'X', '@', '-', '-', '-', '+'],
    ['+', '-', '-', '-', '-', 'X', '-', '-', '+'],
    ['+', 'X', '-', '@', '@', '@', '-', '@', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '-', '@', '@', '-', '-', '-', '+', '+'],
    ['+', '-', '%', '-', '-', '-', '+', '+', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(solution if solution else "No solution found")