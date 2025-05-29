from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, player_pos):
    directions = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    next_states = []
    rows, cols = len(state), len(state[0])
    px, py = player_pos
    
    for dx, dy, move in directions:
        new_px, new_py = px + dx, py + dy
        
        if not is_valid(new_px, new_py, rows, cols) or state[new_px][new_py] == '+':
            continue
            
        new_state = [list(row) for row in state]
        
        # If there's a box at the new position
        if state[new_px][new_py] in ['@', '$']:
            box_x, box_y = new_px + dx, new_py + dy
            
            if not is_valid(box_x, box_y, rows, cols) or state[box_x][box_y] in ['+', '@', '$']:
                continue
                
            # Move the box
            new_state[box_x][box_y] = '$' if state[box_x][box_y] == 'X' else '@'
            new_state[new_px][new_py] = '%' if state[new_px][new_py] == '$' else '-'
        
        # Move the player
        new_state[new_px][new_py] = '%' if state[new_px][new_py] == 'X' else '*'
        new_state[px][py] = 'X' if state[px][py] == '%' else '-'
        
        next_states.append((''.join(''.join(row) for row in new_state), (new_px, new_py), move))
    
    return next_states

def is_solved(state):
    return all(c != '@' for c in ''.join(''.join(row) for row in state))

def solve_sokoban(initial_state):
    rows = len(initial_state)
    cols = len(initial_state[0])
    
    # Find player position
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if initial_state[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    visited = set()
    queue = deque([(initial_state, player_pos, "")])
    visited.add(''.join(''.join(row) for row in initial_state))
    
    while queue:
        current_state, current_pos, path = queue.popleft()
        
        if is_solved(current_state):
            return path
        
        for next_state, next_pos, move in get_next_states(current_state, current_pos):
            if next_state not in visited:
                visited.add(next_state)
                queue.append(([list(next_state[i:i+cols]) for i in range(0, len(next_state), cols)], next_pos, path + move))
    
    return None

# Initial state
initial_state = [
    ['+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '-', '+'],
    ['+', '-', '*', '-', '@', '-', '+'],
    ['+', '-', '-', '@', '$', '+', '+'],
    ['+', '-', 'X', '-', '-', '-', '+'],
    ['+', '+', '$', 'X', '$', '+', '+'],
    ['+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_state)
print(solution)