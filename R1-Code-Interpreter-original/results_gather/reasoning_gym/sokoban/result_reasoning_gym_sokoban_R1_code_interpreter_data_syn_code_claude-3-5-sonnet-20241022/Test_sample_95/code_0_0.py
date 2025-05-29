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
        
        # If moving to empty space or goal
        if state[new_x][new_y] in ['-', 'X']:
            # Update player position
            new_state[player_pos[0]][player_pos[1]] = '-' if state[player_pos[0]][player_pos[1]] == '*' else 'X'
            new_state[new_x][new_y] = '*' if state[new_x][new_y] == '-' else '%'
            next_states.append((new_state, move))
            
        # If moving to box
        elif state[new_x][new_y] in ['@', '$']:
            push_x, push_y = new_x + dx, new_y + dy
            if (is_valid(push_x, push_y, rows, cols) and 
                state[push_x][push_y] in ['-', 'X']):
                # Update box position
                new_state[push_x][push_y] = '@' if state[push_x][push_y] == '-' else '$'
                # Update player position
                new_state[player_pos[0]][player_pos[1]] = '-' if state[player_pos[0]][player_pos[1]] == '*' else 'X'
                new_state[new_x][new_y] = '*' if state[new_x][new_y] == '@' else '%'
                next_states.append((new_state, move))
                
    return next_states

def is_solved(state):
    for row in state:
        for cell in row:
            if cell == '@' or cell == 'X':  # If there's a box not on goal or an empty goal
                return False
    return True

def solve_sokoban(initial_state):
    rows = len(initial_state)
    cols = len(initial_state[0])
    
    queue = deque([(initial_state, "")])
    visited = {tuple(map(tuple, initial_state))}
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
            
        for next_state, move in get_next_states(current_state, rows, cols):
            state_tuple = tuple(map(tuple, next_state))
            if state_tuple not in visited:
                visited.add(state_tuple)
                queue.append((next_state, path + move))
    
    return None

# Initial state
initial_state = [
    list("++++++"),
    list("+-$+++"),
    list("+----+"),
    list("++---+"),
    list("++-X$+"),
    list("+X@@X+"),
    list("+++--+"),
    list("+X-@@+"),
    list("++-*-+"),
    list("++++++")
]

solution = solve_sokoban(initial_state)
print(solution)