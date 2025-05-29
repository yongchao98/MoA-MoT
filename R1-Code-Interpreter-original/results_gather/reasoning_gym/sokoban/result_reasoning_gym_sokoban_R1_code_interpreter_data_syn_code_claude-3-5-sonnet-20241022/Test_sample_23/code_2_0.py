from collections import deque

def get_initial_state():
    return [
        ['+', '+', '+', '+', '+', '+'],
        ['+', '+', 'X', '@', '*', '+'],
        ['+', '+', '-', '-', '@', '+'],
        ['+', '$', '-', '-', 'X', '+'],
        ['+', '+', '+', '+', '+', '+'],
        ['+', '+', '+', '+', '$', '+'],
        ['+', '+', '+', '+', '+', '+']
    ]

def is_goal(state):
    box_count = 0
    for row in state:
        for cell in row:
            if cell == '$':
                box_count += 1
    return box_count == 3

def get_player_pos(state):
    for i in range(len(state)):
        for j in range(len(state[0])):
            if state[i][j] == '*':
                return (i, j)
    return None

def make_move(state, y, x, dy, dx):
    new_state = [list(row) for row in state]
    new_y, new_x = y + dy, x + dx
    
    # Moving to empty space or goal
    if new_state[new_y][new_x] in ['-', 'X']:
        new_state[y][x] = '-'
        new_state[new_y][new_x] = '*'
        return new_state
    
    # Pushing a box
    if new_state[new_y][new_x] in ['@', '$']:
        next_y, next_x = new_y + dy, new_x + dx
        if new_state[next_y][next_x] in ['-', 'X']:
            new_state[y][x] = '-'
            if new_state[next_y][next_x] == 'X':
                new_state[next_y][next_x] = '$'
            else:
                new_state[next_y][next_x] = '@'
            new_state[new_y][new_x] = '*'
            return new_state
    
    return None

def get_state_string(state):
    return '\n'.join(''.join(row) for row in state)

def solve():
    initial = get_initial_state()
    queue = deque([(initial, "")])
    visited = set()
    
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    while queue:
        current_state, path = queue.popleft()
        state_str = get_state_string(current_state)
        
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        if is_goal(current_state):
            return path
        
        y, x = get_player_pos(current_state)
        for move, (dy, dx) in directions.items():
            new_state = make_move(current_state, y, x, dy, dx)
            if new_state:
                new_state_str = get_state_string(new_state)
                if new_state_str not in visited:
                    queue.append((new_state, path + move))
    
    return None

def count_boxes_on_goals(state):
    count = 0
    for row in state:
        count += row.count('$')
    return count

# Solve and print the first move that increases the number of boxes on goals
initial_state = get_initial_state()
initial_boxes_on_goals = count_boxes_on_goals(initial_state)

for move, (dy, dx) in {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}.items():
    y, x = get_player_pos(initial_state)
    new_state = make_move(initial_state, y, x, dy, dx)
    if new_state:
        new_boxes_on_goals = count_boxes_on_goals(new_state)
        if new_boxes_on_goals > initial_boxes_on_goals:
            print(move)
            break