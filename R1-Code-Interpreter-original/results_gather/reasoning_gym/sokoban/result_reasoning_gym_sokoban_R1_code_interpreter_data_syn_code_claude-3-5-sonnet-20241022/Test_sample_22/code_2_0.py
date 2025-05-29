from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state):
    rows, cols = len(state), len(state[0])
    next_states = []
    directions = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    
    # Find all players
    players = []
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['@', '*', '%']:
                players.append((i, j))
    
    for player_pos in players:
        px, py = player_pos
        for dx, dy, move in directions:
            new_x, new_y = px + dx, py + dy
            
            if not is_valid(new_x, new_y, rows, cols):
                continue
                
            if state[new_x][new_y] == '+':
                continue
                
            new_state = [list(row) for row in state]
            
            # If moving to empty space or goal
            if state[new_x][new_y] in ['-', 'X']:
                # Update player position
                new_state[new_x][new_y] = '*' if state[new_x][new_y] == 'X' else '@'
                new_state[px][py] = 'X' if state[px][py] in ['*', '%'] else '-'
                next_states.append((new_state, move))
                
            # If moving a box
            elif state[new_x][new_y] in ['$', '@']:
                push_x, push_y = new_x + dx, new_y + dy
                
                if not is_valid(push_x, push_y, rows, cols):
                    continue
                    
                if state[push_x][push_y] in ['-', 'X']:
                    # Update box position
                    new_state[push_x][push_y] = '$' if state[push_x][push_y] == 'X' else '@'
                    # Update player position
                    new_state[new_x][new_y] = '*' if state[new_x][new_y] == '$' else '@'
                    new_state[px][py] = 'X' if state[px][py] in ['*', '%'] else '-'
                    next_states.append((new_state, move))
    
    return next_states

def is_solved(state):
    for row in state:
        for cell in row:
            if cell == 'X':  # If there's an uncovered goal
                return False
    return True

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = {state_to_string(initial_state)}
    
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