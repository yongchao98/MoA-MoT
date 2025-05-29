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
        
        if not is_valid(new_px, new_py, rows, cols):
            continue
            
        if state[new_px][new_py] in ['+']:
            continue
            
        new_state = [list(row) for row in state]
        
        # If moving to empty space or goal
        if state[new_px][new_py] in ['-', 'X']:
            new_state[px][py] = '-' if state[px][py] == '*' else 'X'
            new_state[new_px][new_py] = '*' if state[new_px][new_py] == '-' else '%'
            next_states.append((new_state, (new_px, new_py), move))
            
        # If pushing a box
        elif state[new_px][new_py] in ['@', '$']:
            box_px, box_py = new_px + dx, new_py + dy
            
            if (is_valid(box_px, box_py, rows, cols) and 
                state[box_px][box_py] in ['-', 'X']):
                
                new_state[px][py] = '-' if state[px][py] == '*' else 'X'
                new_state[new_px][new_py] = '*' if state[new_px][new_py] == '@' else '%'
                new_state[box_px][box_py] = '@' if state[box_px][box_py] == '-' else '$'
                next_states.append((new_state, (new_px, new_py), move))
                
    return next_states

def state_to_string(state):
    return '\n'.join(''.join(row) for row in state)

def solve_sokoban(initial_state):
    # Find player position and convert to list of lists
    initial_state = [list(row) for row in initial_state.split('\n')]
    rows, cols = len(initial_state), len(initial_state[0])
    
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
    visited.add(state_to_string(initial_state))
    
    while queue:
        current_state, current_pos, moves = queue.popleft()
        
        # Check if solved (all boxes are on goals)
        is_solved = True
        for row in current_state:
            if '@' in row:  # If there's still a box not on a goal
                is_solved = False
                break
        
        if is_solved:
            return moves
        
        # Generate next possible states
        for new_state, new_pos, move in get_next_states(current_state, current_pos):
            state_str = state_to_string(new_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((new_state, new_pos, moves + move))
    
    return None

# Initial state
initial_state = """++++++"
++---+"
+--$+"
+X@-++"
+--++"
+@-++"
+%$$++"
++++++"
"""

solution = solve_sokoban(initial_state)
print(solution if solution else "No solution found")