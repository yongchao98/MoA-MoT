from collections import deque

def find_player(state):
    # Return the lowest-rightmost player position
    player_positions = []
    for i in range(len(state)):
        for j in range(len(state[i])):
            if state[i][j] in ['@', '*']:
                player_positions.append((i, j))
    return sorted(player_positions, key=lambda x: (-x[0], -x[1]))[0]

def is_valid_move(state, player_pos, dx, dy):
    rows, cols = len(state), len(state[0])
    new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
    
    if not (0 <= new_x < rows and 0 <= new_y < cols):
        return False
    
    if state[new_x][new_y] == '+':
        return False
        
    if state[new_x][new_y] in ['$', '@']:
        box_x, box_y = new_x + dx, new_y + dy
        if not (0 <= box_x < rows and 0 <= box_y < cols):
            return False
        if state[box_x][box_y] in ['+', '$', '@']:
            return False
    
    return True

def make_move(state, move):
    state = [list(row) for row in state]
    player_pos = find_player(state)
    px, py = player_pos
    
    dx, dy = 0, 0
    if move == 'U': dx = -1
    elif move == 'D': dx = 1
    elif move == 'L': dy = -1
    elif move == 'R': dy = 1
    
    new_x, new_y = px + dx, py + dy
    
    # Handle player's current position
    if state[px][py] == '*':
        state[px][py] = 'X'
    else:
        state[px][py] = '-'
    
    # Handle box movement
    if state[new_x][new_y] in ['$', '@']:
        box_x, box_y = new_x + dx, new_y + dy
        if state[box_x][box_y] == 'X':
            state[box_x][box_y] = '$'
        else:
            state[box_x][box_y] = '@'
    
    # Handle player's new position
    if state[new_x][new_y] == 'X':
        state[new_x][new_y] = '*'
    else:
        state[new_x][new_y] = '@'
    
    return [''.join(row) for row in state]

def get_state_key(state):
    return '\n'.join(state)

def is_solved(state):
    box_count = 0
    goal_count = 0
    box_on_goal = 0
    
    for row in state:
        for cell in row:
            if cell in ['$', '@']:
                box_count += 1
            if cell in ['X', '*']:
                goal_count += 1
            if cell == '$':
                box_on_goal += 1
    
    return box_count == goal_count == box_on_goal

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = {get_state_key(initial_state)}
    moves = ['U', 'D', 'L', 'R']
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
        
        player_pos = find_player(current_state)
        
        for move in moves:
            dx, dy = 0, 0
            if move == 'U': dx = -1
            elif move == 'D': dx = 1
            elif move == 'L': dy = -1
            elif move == 'R': dy = 1
            
            if is_valid_move(current_state, player_pos, dx, dy):
                new_state = make_move(current_state, move)
                new_state_key = get_state_key(new_state)
                
                if new_state_key not in visited:
                    visited.add(new_state_key)
                    queue.append((new_state, path + move))
    
    return None

# Initial state
initial_state = [
    "+++++++" ,
    "+++$+++",
    "+++$+++",
    "++++$$+",
    "+-----+",
    "+----X+",
    "+-X-@-+",
    "+---@*+",
    "+--X@-+",
    "+++++++"
]

solution = solve_sokoban(initial_state)
print(solution)