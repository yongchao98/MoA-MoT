from collections import deque

def find_bottom_player(state):
    # Find the bottom-most player position
    for i in range(len(state)-1, -1, -1):
        for j in range(len(state[i])):
            if state[i][j] in ['@', '*']:
                return (i, j)
    return None

def clean_state(state):
    # Keep only the bottom-most player
    player_pos = find_bottom_player(state)
    new_state = []
    for i in range(len(state)):
        row = list(state[i])
        for j in range(len(row)):
            if (i, j) != player_pos and row[j] in ['@', '*']:
                row[j] = '-'
        new_state.append(''.join(row))
    return new_state

def try_move(state, pos, direction):
    rows, cols = len(state), len(state[0])
    x, y = pos
    dx, dy = 0, 0
    
    if direction == 'U': dx = -1
    elif direction == 'D': dx = 1
    elif direction == 'L': dy = -1
    elif direction == 'R': dy = 1
    
    new_x, new_y = x + dx, y + dy
    
    # Check bounds and walls
    if new_x < 0 or new_x >= rows or new_y < 0 or new_y >= cols or state[new_x][new_y] == '+':
        return None, None
    
    new_state = [list(row) for row in state]
    
    # Moving to empty space or goal
    if state[new_x][new_y] in ['-', 'X']:
        new_state[x][y] = '-' if state[x][y] == '@' else 'X'
        new_state[new_x][new_y] = '@'
        return [''.join(row) for row in new_state], (new_x, new_y)
    
    # Moving box
    if state[new_x][new_y] in ['$', '@']:
        box_x, box_y = new_x + dx, new_y + dy
        if (box_x < 0 or box_x >= rows or box_y < 0 or box_y >= cols or 
            state[box_x][box_y] not in ['-', 'X']):
            return None, None
        
        new_state[box_x][box_y] = '$' if state[box_x][box_y] == 'X' else '@'
        new_state[new_x][new_y] = '@'
        new_state[x][y] = '-' if state[x][y] == '@' else 'X'
        return [''.join(row) for row in new_state], (new_x, new_y)
    
    return None, None

def is_goal_state(state):
    goals = sum(row.count('X') for row in state)
    boxes_on_goals = sum(row.count('$') for row in state)
    return goals == 0 and boxes_on_goals == 3

def solve(initial_state):
    initial_state = clean_state(initial_state)
    start_pos = find_bottom_player(initial_state)
    
    queue = deque([(initial_state, start_pos, "")])
    seen = {''.join(initial_state)}
    
    while queue:
        current_state, pos, path = queue.popleft()
        
        if is_goal_state(current_state):
            return path
        
        for move in ['U', 'D', 'L', 'R']:
            new_state, new_pos = try_move(current_state, pos, move)
            if new_state:
                state_key = ''.join(new_state)
                if state_key not in seen:
                    seen.add(state_key)
                    queue.append((new_state, new_pos, path + move))
    
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

solution = solve(initial_state)
print(solution)