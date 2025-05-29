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
        new_state[new_x][new_y] = '%' if state[new_x][new_y] == 'X' else '*'
        
        next_states.append((''.join(''.join(row) for row in new_state), (new_x, new_y), move))
    
    return next_states

def is_solved(state):
    return all(c != '@' for c in ''.join(''.join(row) for row in state))

def solve_sokoban(initial_state):
    initial_state = [list(row) for row in initial_state.strip().split('\n')]
    
    # Find player position
    player_pos = None
    for i in range(len(initial_state)):
        for j in range(len(initial_state[0])):
            if initial_state[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    visited = set()
    queue = deque([(initial_state, player_pos, "")])
    
    while queue:
        current_state, current_pos, moves = queue.popleft()
        state_str = ''.join(''.join(row) for row in current_state)
        
        if is_solved(current_state):
            return moves
            
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        for next_state, next_pos, move in get_next_states(current_state, current_pos):
            if next_state not in visited:
                queue.append(([list(next_state[i:i+len(current_state[0])]) 
                             for i in range(0, len(next_state), len(current_state[0]))],
                            next_pos, moves + move))
    
    return None

# Initial puzzle state
puzzle = """++++++
+%@X++
+@@-X+
+-@-++
+----+
+--X-+
++++++"""

solution = solve_sokoban(puzzle)
print(solution)