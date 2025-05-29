from collections import deque

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
            
        new_state = [row[:] for row in state]
        
        if state[new_x][new_y] in ['@', '$']:  # Box
            box_x, box_y = new_x + dx, new_y + dy
            if not is_valid(box_x, box_y, rows, cols) or state[box_x][box_y] in ['+', '@', '$']:
                continue
            
            new_state[box_x][box_y] = '$' if (box_x, box_y) in goals else '@'
            new_state[new_x][new_y] = '*'
            new_state[px][py] = '-'
            
        else:  # Empty space or goal
            new_state[new_x][new_y] = '*'
            new_state[px][py] = '-'
            
        next_states.append((new_state, (new_x, new_y), move))
    return next_states

def is_goal_state(state):
    box_count = 0
    goal_count = 0
    boxes_on_goals = 0
    for i in range(len(state)):
        for j in range(len(state[0])):
            if state[i][j] in ['@', '$']:
                box_count += 1
            if state[i][j] in ['X', '$']:
                goal_count += 1
            if state[i][j] == '$':
                boxes_on_goals += 1
    return boxes_on_goals == goal_count == box_count

def solve_sokoban(initial_state):
    rows, cols = len(initial_state), len(initial_state[0])
    
    # Find player position and goals
    player_pos = None
    global goals
    goals = set()
    for i in range(rows):
        for j in range(cols):
            if initial_state[i][j] == '*':
                player_pos = (i, j)
            if initial_state[i][j] == 'X':
                goals.add((i, j))
    
    queue = deque([(initial_state, player_pos, "")])
    visited = {str(initial_state)}
    
    while queue:
        current_state, current_pos, path = queue.popleft()
        
        if is_goal_state(current_state):
            return path
        
        for next_state, next_pos, move in get_next_states(current_state, current_pos):
            state_str = str(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, next_pos, path + move))
    
    return None

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