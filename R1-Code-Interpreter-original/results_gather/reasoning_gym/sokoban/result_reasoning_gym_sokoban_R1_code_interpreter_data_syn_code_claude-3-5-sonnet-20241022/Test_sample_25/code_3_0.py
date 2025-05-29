from collections import deque

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state):
    rows = len(state)
    cols = len(state[0])
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    next_states = []
    
    # Find player position
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if state[i][j] == '*':  # Player
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    if not player_pos:
        return []
        
    for dx, dy, move in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [row[:] for row in state]
        
        # Moving to empty space or goal
        if state[new_x][new_y] in ['-', 'X']:
            new_state[player_pos[0]][player_pos[1]] = '-'  # Previous player position
            new_state[new_x][new_y] = '*'  # New player position
            next_states.append((new_state, move))
            
        # Moving a box
        elif state[new_x][new_y] in ['$', '@']:
            box_new_x, box_new_y = new_x + dx, new_y + dy
            if (is_valid(box_new_x, box_new_y, rows, cols) and 
                state[box_new_x][box_new_y] in ['-', 'X']):
                new_state[player_pos[0]][player_pos[1]] = '-'  # Previous player position
                new_state[new_x][new_y] = '*'  # New player position
                new_state[box_new_x][box_new_y] = '$'  # New box position
                next_states.append((new_state, move))
    
    return next_states

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

def is_solved(state):
    box_count = 0
    goal_count = 0
    box_on_goal = 0
    
    for row in state:
        for cell in row:
            if cell in ['$', '@']:
                box_count += 1
            if cell in ['X', '@']:
                goal_count += 1
            if cell == '@':
                box_on_goal += 1
    
    return box_count == goal_count == box_on_goal

def solve_sokoban(state):
    queue = deque([(state, "")])
    visited = {state_to_string(state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
            
        for next_state, move in get_next_states(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + move))
    
    return None

# Initial state
initial_state = [
    list("++++++"),
    list("+-$+$+"),
    list("+--+++"),
    list("+-@-++"),
    list("+-XX+"),
    list("+-@X@+"),
    list("+--*+"),
    list("++++++")
]

solution = solve_sokoban(initial_state)
print(solution if solution else "No solution found")