from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, rows, cols):
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    next_states = []
    
    # Find player position
    player_pos = None
    boxes = set()
    goals = set()
    
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['*', '@']:
                player_pos = (i, j)
            if state[i][j] in ['$', '@']:
                boxes.add((i, j))
            if state[i][j] in ['X', '@']:
                goals.add((i, j))
    
    for dx, dy, move in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [list(row) for row in state]
        
        # If moving to empty space or goal
        if state[new_x][new_y] in ['-', 'X']:
            new_state[player_pos[0]][player_pos[1]] = '-' if state[player_pos[0]][player_pos[1]] not in ['X', '@'] else 'X'
            new_state[new_x][new_y] = '*' if state[new_x][new_y] == '-' else '%'
            next_states.append((new_state, move))
            
        # If moving a box
        elif state[new_x][new_y] in ['$', '@']:
            box_new_x, box_new_y = new_x + dx, new_y + dy
            if (is_valid(box_new_x, box_new_y, rows, cols) and 
                state[box_new_x][box_new_y] in ['-', 'X']):
                new_state[player_pos[0]][player_pos[1]] = '-' if state[player_pos[0]][player_pos[1]] not in ['X', '@'] else 'X'
                new_state[new_x][new_y] = '*' if state[new_x][new_y] == '$' else '%'
                new_state[box_new_x][box_new_y] = '$' if state[box_new_x][box_new_y] == '-' else '@'
                next_states.append((new_state, move))
    
    return next_states

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

def solve_sokoban(initial_state):
    rows = len(initial_state)
    cols = len(initial_state[0])
    
    queue = deque([(initial_state, "")])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if all boxes are on goals
        all_goals_covered = True
        for i in range(rows):
            for j in range(cols):
                if current_state[i][j] == 'X':
                    all_goals_covered = False
                    break
            if not all_goals_covered:
                break
                
        if all_goals_covered:
            return path
            
        for next_state, move in get_next_states(current_state, rows, cols):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + move))
    
    return None

# Initial state
initial_state = [
    list("++++++"),
    list("+-$+$+"),
    list("+--+++"),
    list("+-@-++"),
    list("+-XX+"),
    list("+-@X@+"),
    list("+--*+"),
    list("++++++")
]

solution = solve_sokoban(initial_state)
print(solution if solution else "No solution found")