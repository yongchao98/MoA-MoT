from collections import deque
import copy

def is_valid_pos(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, rows, cols):
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    
    # Find player position
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    for direction, dx, dy in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid_pos(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [list(row) for row in state]
        
        # Check if moving to empty space or goal
        if state[new_x][new_y] in ['-', 'X']:
            # Update player position
            new_state[player_pos[0]][player_pos[1]] = '-' if state[player_pos[0]][player_pos[1]] == '*' else 'X'
            new_state[new_x][new_y] = '*' if state[new_x][new_y] == '-' else '%'
            next_states.append((direction, new_state))
            
        # Check if pushing a box
        elif state[new_x][new_y] in ['@', '$']:
            push_x, push_y = new_x + dx, new_y + dy
            if (is_valid_pos(push_x, push_y, rows, cols) and 
                state[push_x][push_y] in ['-', 'X']):
                # Update box position
                new_state[push_x][push_y] = '$' if state[push_x][push_y] == 'X' else '@'
                # Update player position
                new_state[player_pos[0]][player_pos[1]] = '-' if state[player_pos[0]][player_pos[1]] == '*' else 'X'
                new_state[new_x][new_y] = '*' if state[new_x][new_y] == '@' else '%'
                next_states.append((direction, new_state))
    
    return next_states

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

def solve_sokoban(initial_state):
    rows = len(initial_state)
    cols = len(initial_state[0])
    
    queue = deque([([], initial_state)])
    visited = {state_to_string(initial_state)}
    
    while queue:
        path, current_state = queue.popleft()
        
        # Check if solved (all boxes on goals)
        if not any('@' in row for row in current_state):
            return ''.join(path)
        
        for direction, next_state in get_next_states(current_state, rows, cols):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((path + [direction], next_state))
    
    return None

# Initial puzzle state
puzzle = [
    ['+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '*', '+'],
    ['+', '@', '@', '@', '@', '+'],
    ['+', 'X', '-', '-', '-', '+'],
    ['+', '-', '-', '-', 'X', '+'],
    ['+', 'X', '-', '$', 'X', '+'],
    ['+', '+', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(puzzle)
print(solution)