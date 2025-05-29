from collections import deque

def get_next_states(state):
    rows = len(state)
    cols = len(state[0])
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    next_states = []
    
    # Find player position
    player_pos = None
    for i in range(rows):
        if '*' in state[i]:
            j = state[i].index('*')
            player_pos = (i, j)
            break
    
    if not player_pos:
        return []
    
    for dx, dy, move in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        # Check bounds
        if new_x < 0 or new_x >= rows or new_y < 0 or new_y >= len(state[new_x]):
            continue
            
        # Check wall
        if state[new_x][new_y] == '+':
            continue
            
        # Create new state
        new_state = []
        for row in state:
            new_state.append(list(row))
            
        # Move to empty space or goal
        if new_state[new_x][new_y] in ['-', 'X']:
            new_state[player_pos[0]][player_pos[1]] = '-'
            new_state[new_x][new_y] = '*'
            next_states.append((new_state, move))
            
        # Push box
        elif new_state[new_x][new_y] in ['$', '@']:
            box_x, box_y = new_x + dx, new_y + dy
            if (0 <= box_x < rows and 
                0 <= box_y < len(state[box_x]) and 
                new_state[box_x][box_y] in ['-', 'X']):
                new_state[player_pos[0]][player_pos[1]] = '-'
                new_state[new_x][new_y] = '*'
                new_state[box_x][box_y] = '$'
                next_states.append((new_state, move))
    
    return next_states

def state_to_string(state):
    return '\n'.join(''.join(row) for row in state)

def is_goal(state):
    goals = 0
    boxes_on_goals = 0
    for row in state:
        for cell in row:
            if cell == 'X':
                goals += 1
            elif cell == '@':
                boxes_on_goals += 1
    return goals == 0 and boxes_on_goals > 0

def solve_sokoban(state):
    queue = deque([(state, "")])
    visited = {state_to_string(state)}
    
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