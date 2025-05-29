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
        
        if state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            
            if not is_valid(box_x, box_y, rows, cols) or state[box_x][box_y] in ['+', '@', '$']:
                continue
                
            # Move box
            new_state[box_x][box_y] = '$' if state[box_x][box_y] == 'X' else '@'
            new_state[new_x][new_y] = '%' if state[new_x][new_y] == '$' else '*'
            new_state[player_pos[0]][player_pos[1]] = 'X' if state[player_pos[0]][player_pos[1]] == '%' else '-'
            
        else:
            new_state[new_x][new_y] = '%' if state[new_x][new_y] == 'X' else '*'
            new_state[player_pos[0]][player_pos[1]] = 'X' if state[player_pos[0]][player_pos[1]] == '%' else '-'
            
        next_states.append((new_state, (new_x, new_y), move))
    
    return next_states

def is_solved(state):
    for row in state:
        for cell in row:
            if cell == '@' or cell == 'X':
                return False
    return True

def solve_sokoban():
    initial_state = [
        ['+', '+', '+', '+', '+', '+', '+'],
        ['+', '+', '-', '-', '-', '-', '+'],
        ['+', '$', '$', '-', '$', '-', '+'],
        ['+', 'X', '-', '-', '-', '-', '+'],
        ['+', '*', '-', '-', '-', '-', '+'],
        ['+', '@', '@', '-', '-', '+', '+'],
        ['+', '-', '@', 'X', 'X', '+', '+'],
        ['+', '-', '-', '-', '-', '-', '+'],
        ['+', '+', '+', '+', '+', '+', '+']
    ]
    
    player_pos = (4, 1)
    visited = set()
    queue = deque([(initial_state, player_pos, "")])
    
    while queue:
        current_state, current_pos, moves = queue.popleft()
        state_hash = str(current_state)
        
        if state_hash in visited:
            continue
            
        visited.add(state_hash)
        
        if is_solved(current_state):
            return moves
            
        for next_state, next_pos, move in get_next_states(current_state, current_pos):
            queue.append((next_state, next_pos, moves + move))
    
    return None

solution = solve_sokoban()
print(f"<<{solution}>>")