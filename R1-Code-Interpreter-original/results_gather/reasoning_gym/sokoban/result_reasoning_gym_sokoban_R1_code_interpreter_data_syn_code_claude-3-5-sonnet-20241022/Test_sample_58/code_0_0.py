from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state):
    rows = len(state)
    cols = len(state[0])
    player_pos = None
    
    # Find player position
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    moves = [(-1,0,'U'), (1,0,'D'), (0,-1,'L'), (0,1,'R')]
    next_states = []
    
    for dx, dy, direction in moves:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [list(row) for row in state]
        
        # If moving to empty space or goal
        if state[new_x][new_y] in ['-', 'X']:
            # Update player position
            new_state[player_pos[0]][player_pos[1]] = '-' if state[player_pos[0]][player_pos[1]] == '*' else 'X'
            new_state[new_x][new_y] = '*' if state[new_x][new_y] == '-' else '%'
            next_states.append((''.join(''.join(row) for row in new_state), direction))
            
        # If moving box
        elif state[new_x][new_y] in ['@', '$']:
            box_new_x, box_new_y = new_x + dx, new_y + dy
            if (is_valid(box_new_x, box_new_y, rows, cols) and 
                state[box_new_x][box_new_y] in ['-', 'X']):
                # Move box
                new_state[box_new_x][box_new_y] = '$' if state[box_new_x][box_new_y] == 'X' else '@'
                # Move player
                new_state[new_x][new_y] = '*' if state[new_x][new_y] == '@' else '%'
                new_state[player_pos[0]][player_pos[1]] = '-' if state[player_pos[0]][player_pos[1]] == '*' else 'X'
                next_states.append((''.join(''.join(row) for row in new_state), direction))
    
    return next_states

def is_solved(state):
    return state.count('$') == 3  # All boxes should be on goals

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = {initial_state}
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
            
        for next_state, direction in get_next_states([list(current_state[i:i+7]) for i in range(0, len(current_state), 7)]):
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + direction))
    
    return None

# Initial state
initial = [
    "+++++++"
    "+*----+"
    "+X@@--+"
    "+X$---+"
    "+-@X--+"
    "+++++++"
]

initial_state = ''.join(initial)
solution = solve_sokoban(initial_state)
print(solution)