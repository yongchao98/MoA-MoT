from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, player_pos):
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    next_states = []
    rows, cols = len(state), len(state[0])
    px, py = player_pos
    
    for dx, dy, move in directions:
        new_x, new_y = px + dx, py + dy
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
            new_state[px][py] = 'X' if state[px][py] == '%' else '-'
            next_states.append((new_state, (new_x, new_y), move))
        else:  # Empty space or goal
            new_state[new_x][new_y] = '%' if state[new_x][new_y] == 'X' else '*'
            new_state[px][py] = 'X' if state[px][py] == '%' else '-'
            next_states.append((new_state, (new_x, new_y), move))
    return next_states

def is_solved(state):
    goals = 0
    boxes_on_goals = 0
    for row in state:
        for cell in row:
            if cell in ['X', '%']:
                goals += 1
            if cell in ['$', '%']:
                boxes_on_goals += 1
    return goals == boxes_on_goals

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

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
    
    queue = deque([(initial_state, player_pos, "")])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, current_pos, path = queue.popleft()
        
        if is_solved(current_state):
            return path
            
        for next_state, next_pos, move in get_next_states(current_state, current_pos):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, next_pos, path + move))
    
    return None

# Initial state
initial_state = [
    ['+', '+', '+', '+', '+', '+'],
    ['+', '+', '-', '*', '-', '+'],
    ['+', '$', '$', '@', '@', '+'],
    ['+', 'X', '-', '-', 'X', '+'],
    ['+', 'X', '@', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '@', 'X', '+'],
    ['+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_state)
print(solution)