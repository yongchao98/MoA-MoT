from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_moves(state):
    rows = len(state)
    cols = len(state[0])
    
    # Find player position
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['@', '*']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    moves = []
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    
    for dx, dy, move in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        if state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if not is_valid(box_x, box_y, rows, cols) or state[box_x][box_y] in ['+', '@', '$']:
                continue
            
        moves.append(move)
    
    return moves

def apply_move(state, move):
    new_state = [list(row) for row in state]
    rows, cols = len(state), len(state[0])
    
    # Find player
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['@', '*']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    dx, dy = 0, 0
    if move == 'R': dy = 1
    elif move == 'L': dy = -1
    elif move == 'D': dx = 1
    elif move == 'U': dx = -1
    
    px, py = player_pos
    new_x, new_y = px + dx, py + dy
    
    # Update player position
    if state[px][py] == '@':
        new_state[px][py] = '-'
    else:  # player on goal
        new_state[px][py] = 'X'
        
    if state[new_x][new_y] in ['$', '@']:
        box_x, box_y = new_x + dx, new_y + dy
        if state[box_x][box_y] == 'X':
            new_state[box_x][box_y] = '$'
        else:
            new_state[box_x][box_y] = '@'
    
    if state[new_x][new_y] == 'X':
        new_state[new_x][new_y] = '*'
    else:
        new_state[new_x][new_y] = '@'
    
    return [''.join(row) for row in new_state]

def is_goal(state):
    for row in state:
        for cell in row:
            if cell == '@':
                return False
    return True

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = {tuple(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal(current_state):
            return path
        
        for move in get_moves(current_state):
            new_state = apply_move(current_state, move)
            state_tuple = tuple(new_state)
            
            if state_tuple not in visited:
                visited.add(state_tuple)
                queue.append((new_state, path + move))
    
    return None

# Initial state
initial_state = [
    "+++++++" ,
    "+++$+++",
    "+++$+++",
    "++++$$+",
    "+-----+",
    "+----X+",
    "+-X-@-+",
    "+---@*+",
    "+--X@-+",
    "+++++++"
]

solution = solve_sokoban(initial_state)
print(solution)