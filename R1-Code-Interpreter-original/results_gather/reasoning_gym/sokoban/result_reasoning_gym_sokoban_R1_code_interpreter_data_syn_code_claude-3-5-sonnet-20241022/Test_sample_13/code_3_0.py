from collections import deque

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_box_positions(state):
    boxes = set()
    for i in range(len(state)):
        for j in range(len(state[0])):
            if state[i][j] in ['@', '$']:
                boxes.add((i, j))
    return boxes

def get_goal_positions(state):
    goals = set()
    for i in range(len(state)):
        for j in range(len(state[0])):
            if state[i][j] in ['X', '$']:
                goals.add((i, j))
    return goals

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def heuristic(state):
    boxes = get_box_positions(state)
    goals = get_goal_positions(state)
    total_distance = 0
    for box in boxes:
        min_dist = float('inf')
        for goal in goals:
            dist = manhattan_distance(box, goal)
            min_dist = min(min_dist, dist)
        total_distance += min_dist
    return total_distance

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
            
            new_state[box_x][box_y] = '$' if state[box_x][box_y] == 'X' else '@'
            new_state[new_x][new_y] = '*'
            new_state[px][py] = '-'
            
        else:  # Empty space or goal
            new_state[new_x][new_y] = '*'
            new_state[px][py] = '-'
            
        next_states.append((new_state, (new_x, new_y), move))
    return next_states

def is_solved(state):
    boxes = get_box_positions(state)
    goals = get_goal_positions(state)
    return boxes == goals

def solve_sokoban(initial_state):
    rows, cols = len(initial_state), len(initial_state[0])
    
    # Find player position
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if initial_state[i][j] == '*':
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    queue = [(heuristic(initial_state), initial_state, player_pos, "")]
    visited = {str(initial_state)}
    
    while queue:
        queue.sort()  # Sort by heuristic value
        _, current_state, current_pos, path = queue.pop(0)
        
        if is_solved(current_state):
            return path
        
        for next_state, next_pos, move in get_next_states(current_state, current_pos):
            state_str = str(next_state)
            if state_str not in visited:
                visited.add(state_str)
                h_value = heuristic(next_state)
                queue.append((h_value, next_state, next_pos, path + move))
    
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