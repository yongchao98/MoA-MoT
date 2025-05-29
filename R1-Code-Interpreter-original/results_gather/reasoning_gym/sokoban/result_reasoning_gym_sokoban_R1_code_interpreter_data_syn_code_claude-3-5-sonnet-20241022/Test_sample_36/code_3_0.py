def get_initial_state(board):
    rows, cols = len(board), len(board[0])
    player = None
    boxes = set()
    goals = set()
    walls = set()
    
    for i in range(rows):
        for j in range(cols):
            if board[i][j] == '+':
                walls.add((i, j))
            elif board[i][j] == '@':
                player = (i, j)
            elif board[i][j] == '$':
                boxes.add((i, j))
            elif board[i][j] == 'X':
                goals.add((i, j))
            elif board[i][j] == '*':
                player = (i, j)
                goals.add((i, j))
            elif board[i][j] == '%':
                player = (i, j)
                goals.add((i, j))
    
    return player, boxes, goals, walls

def is_deadlock(pos, boxes, goals, walls, rows, cols):
    if pos in goals:
        return False
        
    # Corner deadlock
    corners = set()
    for i in range(rows):
        for j in range(cols):
            if (i, j) not in walls:
                continue
            if (i+1, j) in walls and (i, j+1) in walls:
                corners.add((i+1, j+1))
            if (i+1, j) in walls and (i, j-1) in walls:
                corners.add((i+1, j-1))
            if (i-1, j) in walls and (i, j+1) in walls:
                corners.add((i-1, j+1))
            if (i-1, j) in walls and (i, j-1) in walls:
                corners.add((i-1, j-1))
    
    if pos in corners and pos not in goals:
        return True
    
    return False

def get_moves(player, boxes, goals, walls, rows, cols, path, visited):
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    valid_moves = []
    for move, (dx, dy) in directions.items():
        new_player = (player[0] + dx, player[1] + dy)
        
        if new_player in walls:
            continue
            
        if new_player in boxes:
            new_box = (new_player[0] + dx, new_player[1] + dy)
            if new_box not in walls and new_box not in boxes:
                if not is_deadlock(new_box, boxes, goals, walls, rows, cols):
                    new_boxes = set(boxes)
                    new_boxes.remove(new_player)
                    new_boxes.add(new_box)
                    state = (new_player, tuple(sorted(new_boxes)))
                    if state not in visited:
                        valid_moves.append((new_player, new_boxes, path + move))
        else:
            state = (new_player, tuple(sorted(boxes)))
            if state not in visited:
                valid_moves.append((new_player, boxes, path + move))
    
    return valid_moves

def solve_sokoban(board):
    rows, cols = len(board), len(board[0])
    player, boxes, goals, walls = get_initial_state(board)
    
    stack = [(player, boxes, "")]
    visited = {(player, tuple(sorted(boxes)))}
    
    while stack:
        current_player, current_boxes, path = stack.pop()
        
        if current_boxes == goals:
            return path
            
        moves = get_moves(current_player, current_boxes, goals, walls, 
                         rows, cols, path, visited)
        
        for new_player, new_boxes, new_path in moves:
            state = (new_player, tuple(sorted(new_boxes)))
            if state not in visited:
                visited.add(state)
                stack.append((new_player, new_boxes, new_path))
    
    return None

# Initial board
board = [
    ['+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', '+', '+', '-', '$', '+'],
    ['+', '+', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '-', '+'],
    ['+', 'X', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '@', 'X', '-', '+'],
    ['+', '-', '-', '@', '+', '-', '+'],
    ['+', '$', '$', 'X', '-', '@', '+'],
    ['+', '+', '+', '+', '$', '*', '+'],
    ['+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(board)
print(solution)