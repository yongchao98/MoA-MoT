from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, rows, cols):
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    
    # Find player position
    player_x, player_y = None, None
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['%', '*']:
                player_x, player_y = i, j
                break
        if player_x is not None:
            break
    
    for direction, dx, dy in directions:
        new_x, new_y = player_x + dx, player_y + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [list(row) for row in state]
        
        # If moving to empty space or goal
        if state[new_x][new_y] in ['-', 'X']:
            new_state[player_x][player_y] = '-' if state[player_x][player_y] == '%' else 'X'
            new_state[new_x][new_y] = '%' if state[new_x][new_y] == '-' else '*'
            next_states.append((direction, new_state))
            
        # If moving a box
        elif state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if (is_valid(box_x, box_y, rows, cols) and 
                state[box_x][box_y] in ['-', 'X']):
                new_state[player_x][player_y] = '-' if state[player_x][player_y] == '%' else 'X'
                new_state[new_x][new_y] = '%' if state[new_x][new_y] == '@' else '*'
                new_state[box_x][box_y] = '@' if state[box_x][box_y] == '-' else '$'
                next_states.append((direction, new_state))
    
    return next_states

def state_to_string(state):
    return '\n'.join(''.join(row) for row in state)

def is_goal_state(state):
    goals = 0
    for row in state:
        goals += row.count('$')
    return goals == 3

def solve_sokoban(initial_state):
    rows, cols = len(initial_state), len(initial_state[0])
    queue = deque([([], initial_state)])
    visited = {state_to_string(initial_state)}
    
    while queue:
        path, current_state = queue.popleft()
        
        if is_goal_state(current_state):
            return ''.join(path)
            
        for direction, next_state in get_next_states(current_state, rows, cols):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((path + [direction], next_state))
    
    return None

# Initial state
initial_state = [
    ['+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '%', '-', '-', '+'],
    ['+', '-', 'X', '@', '@', '-', '+'],
    ['+', '-', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '-', '@', 'X', '+'],
    ['+', '+', '+', '-', '@', '-', '+'],
    ['+', '+', 'X', '$', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_state)
print(solution)