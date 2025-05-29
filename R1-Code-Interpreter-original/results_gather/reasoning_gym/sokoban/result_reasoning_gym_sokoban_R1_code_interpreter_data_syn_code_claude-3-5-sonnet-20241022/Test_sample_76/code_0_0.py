from collections import deque
import copy

def is_valid(x, y, width, height):
    return 0 <= x < width and 0 <= y < height

def get_moves(state):
    moves = []
    player_pos = None
    width = len(state[0])
    height = len(state)
    
    # Find player position
    for y in range(height):
        for x in range(width):
            if state[y][x] in ['*', '%']:
                player_pos = (x, y)
                break
        if player_pos:
            break
    
    directions = [('L', -1, 0), ('R', 1, 0), ('U', 0, -1), ('D', 0, 1)]
    
    for direction, dx, dy in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, width, height) or state[new_y][new_x] == '+':
            continue
            
        if state[new_y][new_x] in ['@', '$']:
            box_new_x, box_new_y = new_x + dx, new_y + dy
            if (not is_valid(box_new_x, box_new_y, width, height) or 
                state[box_new_y][box_new_x] in ['+', '@', '$']):
                continue
            moves.append((direction, (new_x, new_y), (box_new_x, box_new_y)))
        else:
            moves.append((direction, (new_x, new_y), None))
    
    return moves

def apply_move(state, move):
    new_state = [list(row) for row in state]
    direction, new_player_pos, new_box_pos = move
    old_player_pos = None
    
    # Find current player position
    for y in range(len(state)):
        for x in range(len(state[0])):
            if state[y][x] in ['*', '%']:
                old_player_pos = (x, y)
                break
        if old_player_pos:
            break
    
    # Update player position
    if state[old_player_pos[1]][old_player_pos[0]] == '%':
        new_state[old_player_pos[1]][old_player_pos[0]] = 'X'
    else:
        new_state[old_player_pos[1]][old_player_pos[0]] = '-'
    
    if state[new_player_pos[1]][new_player_pos[0]] in ['X', '$']:
        new_state[new_player_pos[1]][new_player_pos[0]] = '%'
    else:
        new_state[new_player_pos[1]][new_player_pos[0]] = '*'
    
    # Update box position if applicable
    if new_box_pos:
        if state[new_player_pos[1]][new_player_pos[0]] == '$':
            new_state[new_player_pos[1]][new_player_pos[0]] = '%'
        else:
            new_state[new_player_pos[1]][new_player_pos[0]] = '*'
            
        if state[new_box_pos[1]][new_box_pos[0]] == 'X':
            new_state[new_box_pos[1]][new_box_pos[0]] = '$'
        else:
            new_state[new_box_pos[1]][new_box_pos[0]] = '@'
    
    return [''.join(row) for row in new_state]

def is_solved(state):
    goals = 0
    boxes_on_goals = 0
    for row in state:
        for cell in row:
            if cell in ['X', '%']:
                goals += 1
            if cell in ['$', '%']:
                boxes_on_goals += 1
    return goals == boxes_on_goals

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    seen = {tuple(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
        
        for move in get_moves(current_state):
            new_state = apply_move(current_state, move)
            new_state_tuple = tuple(new_state)
            
            if new_state_tuple not in seen:
                seen.add(new_state_tuple)
                queue.append((new_state, path + move[0]))

# Initial state
initial_state = [
    "+++++++",
    "++X@*+",
    "++XX@-+",
    "++@-+",
    "++$+-+",
    "+++++++"
]

solution = solve_sokoban(initial_state)
print(solution)