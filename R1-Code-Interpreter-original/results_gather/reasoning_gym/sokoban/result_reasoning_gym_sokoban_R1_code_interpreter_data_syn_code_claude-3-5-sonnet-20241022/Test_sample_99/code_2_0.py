from collections import deque

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, rows, cols, active_player):
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    
    # Find the active player position
    px = py = None
    for i in range(rows):
        for j in range(cols):
            if (i, j) == active_player:
                px, py = i, j
                break
        if px is not None:
            break
    
    for move, dx, dy in directions:
        new_x, new_y = px + dx, py + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [row[:] for row in state]
        
        # Moving to empty space or goal
        if state[new_x][new_y] in ['-', 'X']:
            new_state[new_x][new_y] = '@'
            new_state[px][py] = '-'
            next_states.append((move, new_state, (new_x, new_y)))
            
        # Moving a box
        elif state[new_x][new_y] in ['$', '@']:
            box_x, box_y = new_x + dx, new_y + dy
            if (is_valid(box_x, box_y, rows, cols) and 
                state[box_x][box_y] in ['-', 'X']):
                new_state[box_x][box_y] = '$' if state[box_x][box_y] == 'X' else '@'
                new_state[new_x][new_y] = '@'
                new_state[px][py] = '-'
                next_states.append((move, new_state, (new_x, new_y)))
    
    return next_states

def is_solved(state):
    boxes_on_goals = sum(row.count('$') for row in state)
    return boxes_on_goals == 2

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

def solve_sokoban(initial_state):
    # Start with the leftmost player
    start_pos = None
    for i in range(len(initial_state)):
        for j in range(len(initial_state[0])):
            if initial_state[i][j] in ['@', '%']:
                start_pos = (i, j)
                break
        if start_pos:
            break
    
    queue = deque([([], initial_state, start_pos)])
    seen = set()
    max_moves = 15
    
    while queue:
        moves, current_state, player_pos = queue.popleft()
        
        if len(moves) > max_moves:
            continue
        
        state_str = state_to_string(current_state)
        if state_str in seen:
            continue
        seen.add(state_str)
        
        if is_solved(current_state):
            return ''.join(moves)
        
        for move, next_state, new_player_pos in get_next_states(current_state, 
                                                               len(current_state), 
                                                               len(current_state[0]),
                                                               player_pos):
            queue.append((moves + [move], next_state, new_player_pos))
    
    return None

# Initial state
initial_state = [
    ['+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', '+', '$', 'X', '@', '%', '+'],
    ['+', '+', '$', '+', '@', '@', '@', '+'],
    ['+', '+', '+', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '$', '-', '-', '+', '+'],
    ['+', '+', '+', '+', '+', 'X', '+', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_state)
print(solution)