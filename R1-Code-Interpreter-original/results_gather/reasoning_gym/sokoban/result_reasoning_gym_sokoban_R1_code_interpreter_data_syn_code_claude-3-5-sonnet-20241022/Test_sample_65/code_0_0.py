from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, rows, cols):
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
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
    
    for dx, dy, move in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [list(row) for row in state]
        
        # If next position is a box
        if state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            
            if not is_valid(box_x, box_y, rows, cols) or \
               state[box_x][box_y] in ['+', '@', '$']:
                continue
                
            # Move box
            if state[box_x][box_y] == 'X':
                new_state[box_x][box_y] = '$'
            else:
                new_state[box_x][box_y] = '@'
            
            # Update previous box position
            if state[new_x][new_y] == '$':
                new_state[new_x][new_y] = 'X'
            else:
                new_state[new_x][new_y] = '-'
        
        # Update player position
        if state[player_pos[0]][player_pos[1]] == '%':
            new_state[player_pos[0]][player_pos[1]] = 'X'
        else:
            new_state[player_pos[0]][player_pos[1]] = '-'
            
        if state[new_x][new_y] == 'X' or \
           (state[new_x][new_y] == '$' and state[player_pos[0]][player_pos[1]] == '%'):
            new_state[new_x][new_y] = '%'
        else:
            new_state[new_x][new_y] = '*'
            
        next_states.append((new_state, move))
    
    return next_states

def state_to_string(state):
    return '\n'.join(''.join(row) for row in state)

def solve_sokoban(initial_state):
    rows = len(initial_state)
    cols = len(initial_state[0])
    
    queue = deque([(initial_state, "")])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if solved
        is_solved = True
        for i in range(rows):
            for j in range(cols):
                if current_state[i][j] == 'X':
                    is_solved = False
                    break
            if not is_solved:
                break
        
        if is_solved:
            return path
        
        for next_state, move in get_next_states(current_state, rows, cols):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + move))
    
    return None

# Initial state
initial_state = [
    ['+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', '+', 'X', '$', '+', '+'],
    ['+', '+', '-', '-', '-', '-', '+'],
    ['+', 'X', '-', '@', '@', '-', '+'],
    ['+', 'X', '-', '-', '-', '-', '+'],
    ['+', 'X', '-', '@', '-', '-', '+'],
    ['+', '@', '-', '@', '-', '-', '+'],
    ['+', 'X', '-', '-', '*', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_state)
print(f"<<<{solution}>>>")