from collections import deque

def get_next_states(state):
    rows, cols = len(state), len(state[0])
    moves = []
    
    # Find players and boxes
    players = []
    boxes = set()
    goals = set()
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['*', '%']:
                players.append((i, j))
            elif state[i][j] == '@':
                if len(players) < 2:  # Second @ is a player
                    players.append((i, j))
                else:  # Other @ are boxes
                    boxes.add((i, j))
            elif state[i][j] == '$':
                boxes.add((i, j))
            elif state[i][j] == 'X':
                goals.add((i, j))
    
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    
    for px, py in players:
        for dx, dy, direction in directions:
            new_x, new_y = px + dx, py + dy
            
            if not (0 <= new_x < rows and 0 <= new_y < cols) or state[new_x][new_y] == '+':
                continue
                
            new_state = [row[:] for row in state]
            current_cell = new_state[px][py]
            target_cell = new_state[new_x][new_y]
            
            # Moving to empty space
            if target_cell in ['-', 'X']:
                new_state[new_x][new_y] = '*'
                new_state[px][py] = 'X' if (px, py) in goals else '-'
                moves.append((new_state, direction))
            
            # Moving a box
            elif target_cell in ['@', '$']:
                push_x, push_y = new_x + dx, new_y + dy
                if (0 <= push_x < rows and 0 <= push_y < cols and 
                    new_state[push_x][push_y] in ['-', 'X']):
                    # Move box
                    new_state[push_x][push_y] = '$' if new_state[push_x][push_y] == 'X' else '@'
                    # Move player
                    new_state[new_x][new_y] = '*'
                    # Update old player position
                    new_state[px][py] = 'X' if (px, py) in goals else '-'
                    moves.append((new_state, direction))
    
    return moves

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

def is_goal(state):
    goals = 0
    boxes_on_goals = 0
    for i, row in enumerate(state):
        for j, cell in enumerate(row):
            if cell == 'X' or cell == '$':
                goals += 1
            if cell == '$':
                boxes_on_goals += 1
    return goals == 3 and boxes_on_goals == 3

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = {state_to_string(initial_state)}
    max_moves = 30  # Limit search depth
    
    while queue:
        current_state, path = queue.popleft()
        
        if len(path) > max_moves:
            continue
            
        if is_goal(current_state):
            return path
            
        for next_state, move in get_next_states(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + move))
                
    return None

# Initial state
initial_state = [
    ['+', '+', '+', '+', '+', '+'],
    ['+', '+', 'X', '@', '*', '+'],
    ['+', '+', '+', '$', '@', '+'],
    ['+', '+', '+', '+', '-', '+'],
    ['+', '+', '$', '+', '-', '+'],
    ['+', '+', '+', '+', '-', '+'],
    ['+', '+', '$', '+', 'X', '+'],
    ['+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_state)
print(solution)