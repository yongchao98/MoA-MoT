from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, player_pos):
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    next_states = []
    rows, cols = len(state), len(state[0])
    
    for dx, dy, move in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [list(row) for row in state]
        new_player_pos = (new_x, new_y)
        
        if state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if not is_valid(box_x, box_y, rows, cols) or state[box_x][box_y] in ['+', '@', '$']:
                continue
                
            # Move box
            if state[new_x][new_y] == '@':
                new_state[new_x][new_y] = '-'
            else:  # '$'
                new_state[new_x][new_y] = 'X'
                
            if state[box_x][box_y] == 'X':
                new_state[box_x][box_y] = '$'
            else:  # '-'
                new_state[box_x][box_y] = '@'
        
        # Move player
        if state[player_pos[0]][player_pos[1]] == '*':
            new_state[player_pos[0]][player_pos[1]] = '-'
        else:  # '%'
            new_state[player_pos[0]][player_pos[1]] = 'X'
            
        if state[new_x][new_y] in ['X', '$']:
            new_state[new_x][new_y] = '%'
        else:
            new_state[new_x][new_y] = '*'
            
        next_states.append((new_state, new_player_pos, move))
    return next_states

def is_solved(state):
    for row in state:
        for cell in row:
            if cell in ['@', 'X']:  # If there's a box not on goal or goal without box
                return False
    return True

def state_to_string(state):
    return '\n'.join(''.join(row) for row in state)

def solve_sokoban(initial_state):
    # Find player position
    player_pos = None
    for i in range(len(initial_state)):
        for j in range(len(initial_state[0])):
            if initial_state[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    queue = deque([(initial_state, player_pos, "")])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, current_pos, moves = queue.popleft()
        
        if is_solved(current_state):
            return moves
            
        for new_state, new_pos, move in get_next_states(current_state, current_pos):
            state_str = state_to_string(new_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((new_state, new_pos, moves + move))
    
    return None

# Initial state
initial_state = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', 'X', '-', '@', 'X', '-', '@', '-', '+'],
    ['+', '@', 'X', '@', '-', '@', '@', 'X', '+'],
    ['+', '-', '-', '-', '-', 'X', '-', '-', '+'],
    ['+', '-', '*', '@', 'X', '-', '-', '-', '+'],
    ['+', '-', '-', '-', 'X', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+']]

solution = solve_sokoban(initial_state)
print(solution if solution else "No solution found")