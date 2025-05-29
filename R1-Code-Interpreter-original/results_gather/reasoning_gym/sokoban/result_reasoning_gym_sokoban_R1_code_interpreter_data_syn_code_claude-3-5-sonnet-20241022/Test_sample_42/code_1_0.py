def get_initial_state():
    puzzle = [
        ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'],
        ['+', '-', '-', 'X', '-', '-', '-', '@', '*', '+'],
        ['+', '-', '-', '-', '@', 'X', '@', '@', '@', '+'],
        ['+', '-', '-', '-', 'X', '+', 'X', '$', '-', '+'],
        ['+', '+', '-', '-', '+', '$', '+', '+', 'X', '+'],
        ['+', 'X', '-', '-', '+', '+', '+', '+', '$', '+'],
        ['+', '-', '@', '-', '+', '+', '+', '+', '+', '+'],
        ['+', '-', '-', '-', '-', '+', '+', '+', '+', '+'],
        ['+', '-', '-', '-', '-', '+', '+', '+', '$', '+'],
        ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+']]
    return puzzle

def get_player_pos(state):
    for i in range(len(state)):
        for j in range(len(state[i])):
            if state[i][j] in ['@', '*']:
                return (i, j)
    return None

def is_box(char):
    return char in ['@', '$']

def is_goal(char):
    return char in ['X', '*', '$']

def can_move(state, r, c):
    return state[r][c] not in ['+']

def make_move(state, move, player_pos):
    new_state = [row[:] for row in state]
    r, c = player_pos
    
    if move == 'U': dr, dc = -1, 0
    elif move == 'D': dr, dc = 1, 0
    elif move == 'L': dr, dc = 0, -1
    else: dr, dc = 0, 1
    
    new_r, new_c = r + dr, c + dc
    
    if not (0 <= new_r < len(state) and 0 <= new_c < len(state[0])):
        return None
        
    curr_char = state[r][c]
    next_char = state[new_r][new_c]
    
    # Moving to empty space or goal
    if next_char in ['-', 'X']:
        new_state[r][c] = 'X' if curr_char == '*' else '-'
        new_state[new_r][new_c] = '*' if next_char == 'X' else '@'
        return new_state
        
    # Pushing a box
    elif is_box(next_char):
        push_r, push_c = new_r + dr, new_c + dc
        if (0 <= push_r < len(state) and 0 <= push_c < len(state[0])):
            push_char = state[push_r][push_c]
            if push_char in ['-', 'X']:
                new_state[r][c] = 'X' if curr_char == '*' else '-'
                new_state[new_r][new_c] = '*' if next_char == '$' else '@'
                new_state[push_r][push_c] = '$' if push_char == 'X' else '@'
                return new_state
    
    return None

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

def is_goal_state(state):
    boxes = 0
    boxes_on_goals = 0
    for row in state:
        for cell in row:
            if cell in ['@', '$', '*']:
                boxes += 1
            if cell in ['$', '*']:
                boxes_on_goals += 1
    return boxes == boxes_on_goals

def solve_iddfs(state, depth_limit, path="", visited=None):
    if visited is None:
        visited = set()
    
    if len(path) > depth_limit:
        return None
        
    if is_goal_state(state):
        return path
        
    state_str = state_to_string(state)
    if state_str in visited:
        return None
        
    visited.add(state_str)
    player_pos = get_player_pos(state)
    
    for move in ['U', 'D', 'L', 'R']:
        new_state = make_move(state, move, player_pos)
        if new_state:
            result = solve_iddfs(new_state, depth_limit, path + move, visited)
            if result:
                return result
                
    visited.remove(state_str)
    return None

def solve_sokoban():
    initial_state = get_initial_state()
    
    # Try increasingly deeper searches
    for depth in range(1, 40):
        result = solve_iddfs(initial_state, depth)
        if result:
            return result
    
    return None

solution = solve_sokoban()
print(solution if solution else "No solution found")