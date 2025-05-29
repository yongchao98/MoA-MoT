from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, player_pos):
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    rows, cols = len(state), len(state[0])
    
    for direction, dx, dy in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [list(row) for row in state]
        if state[new_x][new_y] in ['@', '$']:  # Box
            box_x, box_y = new_x + dx, new_y + dy
            if not is_valid(box_x, box_y, rows, cols) or state[box_x][box_y] in ['+', '@', '$']:
                continue
                
            # Move box
            new_state[box_x][box_y] = '$' if state[box_x][box_y] == 'X' else '@'
            new_state[new_x][new_y] = '%' if state[new_x][new_y] == '$' else '*'
            new_state[player_pos[0]][player_pos[1]] = 'X' if state[player_pos[0]][player_pos[1]] == '%' else '-'
            
        else:  # Empty space or goal
            new_state[new_x][new_y] = '%' if state[new_x][new_y] == 'X' else '*'
            new_state[player_pos[0]][player_pos[1]] = 'X' if state[player_pos[0]][player_pos[1]] == '%' else '-'
            
        next_states.append((direction, new_state, (new_x, new_y)))
    
    return next_states

def is_solved(state):
    for row in state:
        for cell in row:
            if cell in ['@', 'X']:  # If there's a box not on goal or an unfilled goal
                return False
    return True

def solve_sokoban(initial_state):
    # Find player position
    player_pos = None
    for i in range(len(initial_state)):
        for j in range(len(initial_state[0])):
            if initial_state[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    queue = deque([([], initial_state, player_pos)])
    visited = set()
    
    while queue:
        path, current_state, pos = queue.popleft()
        state_str = ''.join(''.join(row) for row in current_state)
        
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        if is_solved(current_state):
            return ''.join(path)
            
        for direction, new_state, new_pos in get_next_states(current_state, pos):
            queue.append((path + [direction], new_state, new_pos))
    
    return None

# Initial state
initial_state = [
    list("+++++++"),
    list("++--@X++"),
    list("++----+"),
    list("++++$-+"),
    list("+++-$-+"),
    list("+$+---+"),
    list("+*@-X-+"),
    list("+++++++")
]

solution = solve_sokoban(initial_state)
print(solution)