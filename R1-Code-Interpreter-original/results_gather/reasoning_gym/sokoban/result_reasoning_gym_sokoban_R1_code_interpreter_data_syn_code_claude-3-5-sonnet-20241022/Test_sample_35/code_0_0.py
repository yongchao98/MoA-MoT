from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_moves(state):
    moves = []
    player_pos = None
    rows = len(state)
    cols = len(state[0])
    
    # Find player position
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    directions = [(-1,0,'U'), (1,0,'D'), (0,-1,'L'), (0,1,'R')]
    
    for dx, dy, move in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols):
            continue
            
        if state[new_x][new_y] in ['+']:
            continue
            
        # If next position is a box
        if state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if not is_valid(box_x, box_y, rows, cols):
                continue
            if state[box_x][box_y] not in ['-', 'X']:
                continue
            
            # Create new state
            new_state = [list(row) for row in state]
            
            # Move box
            if state[box_x][box_y] == 'X':
                new_state[box_x][box_y] = '$'
            else:
                new_state[box_x][box_y] = '@'
                
            # Move player
            if state[new_x][new_y] == '$':
                new_state[new_x][new_y] = '%'
            else:
                new_state[new_x][new_y] = '*'
                
            # Clear old player position
            if state[player_pos[0]][player_pos[1]] == '%':
                new_state[player_pos[0]][player_pos[1]] = 'X'
            else:
                new_state[player_pos[0]][player_pos[1]] = '-'
                
            moves.append((''.join(''.join(row) for row in new_state), move))
        else:
            # Simple player move
            new_state = [list(row) for row in state]
            
            # Move player
            if state[new_x][new_y] == 'X':
                new_state[new_x][new_y] = '%'
            else:
                new_state[new_x][new_y] = '*'
                
            # Clear old player position
            if state[player_pos[0]][player_pos[1]] == '%':
                new_state[player_pos[0]][player_pos[1]] = 'X'
            else:
                new_state[player_pos[0]][player_pos[1]] = '-'
                
            moves.append((''.join(''.join(row) for row in new_state), move))
    
    return moves

def is_solved(state):
    return '@' not in state and '*' not in state

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = {initial_state}
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
            
        for new_state, move in get_moves(current_state.split('\n')):
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + move))
    
    return None

# Initial state
initial = """++++++
+*$--+
+@@-++
+XX+-+
+++$-+
++++++"""

solution = solve_sokoban(initial)
print(solution)