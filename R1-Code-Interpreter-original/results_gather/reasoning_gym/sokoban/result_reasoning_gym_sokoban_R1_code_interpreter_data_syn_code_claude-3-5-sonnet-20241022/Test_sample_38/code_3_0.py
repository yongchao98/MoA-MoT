from collections import deque

def get_next_states(state):
    rows, cols = len(state), len(state[0])
    moves = []
    
    # Find players positions
    players = []
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['*', '@']:
                players.append((i, j))
    
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    
    for px, py in players:
        for dx, dy, direction in directions:
            new_x, new_y = px + dx, py + dy
            
            # Check if move is within bounds
            if not (0 <= new_x < rows and 0 <= new_y < cols):
                continue
            
            # Create new state
            new_state = [row[:] for row in state]
            
            # Moving to empty space or goal
            if state[new_x][new_y] in ['-', 'X']:
                new_state[new_x][new_y] = '*'
                new_state[px][py] = '-'
                moves.append((new_state, direction))
            
            # Moving a box
            elif state[new_x][new_y] in ['@', '$']:
                push_x, push_y = new_x + dx, new_y + dy
                if (0 <= push_x < rows and 0 <= push_y < cols and 
                    state[push_x][push_y] in ['-', 'X']):
                    new_state[push_x][push_y] = '$' if state[push_x][push_y] == 'X' else '@'
                    new_state[new_x][new_y] = '*'
                    new_state[px][py] = '-'
                    moves.append((new_state, direction))
    
    return moves

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

def is_goal(state):
    box_count = sum(row.count('@') + row.count('$') for row in state)
    boxes_on_goals = sum(row.count('$') for row in state)
    return box_count > 0 and box_count == boxes_on_goals

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
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