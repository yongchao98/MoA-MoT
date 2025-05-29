from collections import deque

def get_initial_state():
    return [
        ['+', '+', '+', '+', '+', '+'],
        ['+', '+', 'X', '@', '*', '+'],
        ['+', '+', '-', '-', '@', '+'],
        ['+', '$', '-', '-', 'X', '+'],
        ['+', '+', '+', '+', '+', '+'],
        ['+', '+', '+', '+', '$', '+'],
        ['+', '+', '+', '+', '+', '+']
    ]

def is_valid_pos(y, x, state):
    return 0 <= y < len(state) and 0 <= x < len(state[0])

def get_player_pos(state):
    for i in range(len(state)):
        for j in range(len(state[0])):
            if state[i][j] == '*':
                return (i, j)
    return None

def can_move(state, y, x, dy, dx):
    new_y, new_x = y + dy, x + dx
    if not is_valid_pos(new_y, new_x, state):
        return False
    if state[new_y][new_x] == '+':
        return False
    if state[new_y][new_x] in ['@', '$']:
        next_y, next_x = new_y + dy, new_x + dx
        if not is_valid_pos(next_y, next_x, state):
            return False
        if state[next_y][next_x] not in ['-', 'X']:
            return False
    return True

def make_move(state, move):
    state = [list(row) for row in state]
    y, x = get_player_pos(state)
    
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    dy, dx = directions[move]
    new_y, new_x = y + dy, x + dx
    
    # Move player
    state[y][x] = '-'
    
    # If pushing box
    if state[new_y][new_x] in ['@', '$']:
        next_y, next_x = new_y + dy, new_x + dx
        # Move box
        if state[next_y][next_x] == 'X':
            state[next_y][next_x] = '$'
        else:
            state[next_y][next_x] = '@'
        state[new_y][new_x] = '*'
    else:
        state[new_y][new_x] = '*'
    
    return [''.join(row) for row in state]

def get_state_string(state):
    return ''.join(''.join(row) for row in state)

def solve():
    initial = get_initial_state()
    queue = deque([(initial, "")])
    visited = set()
    
    while queue:
        state, path = queue.popleft()
        state_str = get_state_string(state)
        
        if state_str in visited:
            continue
        
        visited.add(state_str)
        
        y, x = get_player_pos(state)
        for move, (dy, dx) in {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}.items():
            if can_move(state, y, x, dy, dx):
                new_state = make_move(state, move)
                new_path = path + move
                
                # Check if we've pushed a box to a goal
                if any(row.count('$') > row.count('$') for row, new_row in zip(state, new_state)):
                    return new_path
                
                queue.append((new_state, new_path))
    
    return None

result = solve()
print(result)