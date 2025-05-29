from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, player_pos):
    rows, cols = len(state), len(state[0])
    moves = [(-1,0,'U'), (1,0,'D'), (0,-1,'L'), (0,1,'R')]
    next_states = []
    
    for dx, dy, direction in moves:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [list(row) for row in state]
        
        if state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if not is_valid(box_x, box_y, rows, cols) or state[box_x][box_y] in ['+', '@', '$']:
                continue
            
            # Move box
            if state[box_x][box_y] == 'X':
                new_state[box_x][box_y] = '$'
            else:
                new_state[box_x][box_y] = '@'
            
            if state[new_x][new_y] == '$':
                new_state[new_x][new_y] = 'X'
            else:
                new_state[new_x][new_y] = '-'
        
        # Move player
        if state[player_pos[0]][player_pos[1]] == '%':
            new_state[player_pos[0]][player_pos[1]] = 'X'
        else:
            new_state[player_pos[0]][player_pos[1]] = '-'
            
        if state[new_x][new_y] == 'X':
            new_state[new_x][new_y] = '%'
        else:
            new_state[new_x][new_y] = '*'
            
        next_states.append(([''.join(row) for row in new_state], (new_x, new_y), direction))
    
    return next_states

def is_solved(state):
    return all(row.count('@') == 0 for row in state)

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
    visited.add(tuple(initial_state))
    
    while queue:
        current_state, current_pos, path = queue.popleft()
        
        if is_solved(current_state):
            return path
        
        for next_state, next_pos, direction in get_next_states(current_state, current_pos):
            state_tuple = tuple(next_state)
            if state_tuple not in visited:
                visited.add(state_tuple)
                queue.append((next_state, next_pos, path + direction))
    
    return None

# Initial state
initial_state = [
    "+++++++++",
    "+X-*----+",
    "+--@X-@@+",
    "+---++--+",
    "+++$$X@-+",
    "++++++-X+",
    "+++++++++"
]

solution = solve_sokoban(initial_state)
print(solution)