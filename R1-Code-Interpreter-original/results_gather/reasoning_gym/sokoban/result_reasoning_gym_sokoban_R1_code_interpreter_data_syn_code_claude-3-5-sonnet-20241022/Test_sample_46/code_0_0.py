from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state):
    rows, cols = len(state), len(state[0])
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    directions = [(-1,0,'U'), (1,0,'D'), (0,-1,'L'), (0,1,'R')]
    next_states = []
    
    for dx, dy, move in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [list(row) for row in state]
        
        if state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if not is_valid(box_x, box_y, rows, cols) or state[box_x][box_y] in ['+', '@', '$']:
                continue
                
            # Move box
            new_state[box_x][box_y] = '$' if state[box_x][box_y] == 'X' else '@'
            new_state[new_x][new_y] = '%' if state[new_x][new_y] == '$' else '*'
            
        else:
            new_state[new_x][new_y] = '%' if state[new_x][new_y] == 'X' else '*'
            
        new_state[player_pos[0]][player_pos[1]] = 'X' if state[player_pos[0]][player_pos[1]] == '%' else '-'
        next_states.append((''.join(''.join(row) for row in new_state), move))
    
    return next_states

def is_solved(state):
    return all(c != '@' for c in ''.join(''.join(row) for row in state))

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set(''.join(''.join(row) for row in initial_state))
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
            
        for next_state, move in get_next_states(current_state):
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + move))
    
    return None

# Initial state
initial_state = [
    list("+++++++""),
    list("++*$--+""),
    list("+@-@X-+""),
    list("+-----+""),
    list("+X+-@X+""),
    list("+++++++"")
]

solution = solve_sokoban(initial_state)
print(f"<<<{solution}>>>")