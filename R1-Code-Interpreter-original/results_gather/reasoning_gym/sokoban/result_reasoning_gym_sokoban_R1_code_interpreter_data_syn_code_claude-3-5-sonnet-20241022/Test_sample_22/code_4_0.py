from collections import deque

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state):
    rows, cols = len(state), len(state[0])
    next_states = []
    directions = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    
    # Find the active player (prioritize the one that can make meaningful moves)
    players = []
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['@', '*']:
                players.append((i, j))
    
    for px, py in players:
        for dx, dy, move in directions:
            new_x, new_y = px + dx, py + dy
            
            if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
                continue
                
            new_state = [row[:] for row in state]
            
            # Moving to empty space or goal
            if state[new_x][new_y] in ['-', 'X']:
                new_state[new_x][new_y] = '*' if state[new_x][new_y] == 'X' else '@'
                new_state[px][py] = 'X' if state[px][py] == '*' else '-'
                next_states.append((new_state, move))
            
            # Moving a box
            elif state[new_x][new_y] == '$':
                push_x, push_y = new_x + dx, new_y + dy
                
                if (not is_valid(push_x, push_y, rows, cols) or 
                    state[push_x][push_y] not in ['-', 'X']):
                    continue
                
                new_state[push_x][push_y] = '$'
                new_state[new_x][new_y] = '*' if state[new_x][new_y] == '$' else '@'
                new_state[px][py] = 'X' if state[px][py] == '*' else '-'
                next_states.append((new_state, move))
    
    return next_states

def count_goals_and_boxes(state):
    goals = 0
    boxes_on_goals = 0
    for row in state:
        for cell in row:
            if cell == 'X':
                goals += 1
            elif cell == '$':
                boxes_on_goals += 1
    return goals, boxes_on_goals

def is_solved(state):
    goals, boxes_on_goals = count_goals_and_boxes(state)
    return goals == 0 and boxes_on_goals == 3

def state_to_string(state):
    return '\n'.join(''.join(row) for row in state)

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = {state_to_string(initial_state)}
    max_depth = 30  # Limit search depth to avoid infinite loops
    
    while queue:
        current_state, path = queue.popleft()
        
        if len(path) > max_depth:
            continue
            
        if is_solved(current_state):
            return path
            
        for next_state, move in get_next_states(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + move))
    
    return None

initial_state = [
    ['+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '$', '-', '+'],
    ['+', '-', 'X', '-', '@', '-', '+'],
    ['+', '-', '-', '-', '*', '-', '+'],
    ['+', '-', '$', 'X', 'X', '-', '+'],
    ['+', '+', '-', '-', '$', '-', '+'],
    ['+', '-', '@', '@', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_state)
print(solution)