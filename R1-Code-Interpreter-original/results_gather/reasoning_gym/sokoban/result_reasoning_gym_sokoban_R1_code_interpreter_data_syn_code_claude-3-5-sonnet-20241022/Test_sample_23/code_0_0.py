from collections import deque
import copy

def get_initial_state():
    state = [
        ['+', '+', '+', '+', '+', '+'],
        ['+', '+', 'X', '@', '*', '+'],
        ['+', '+', '-', '-', '@', '+'],
        ['+', '$', '-', '-', 'X', '+'],
        ['+', '+', '+', '+', '+', '+'],
        ['+', '+', '+', '+', '$', '+'],
        ['+', '+', '+', '+', '+', '+']
    ]
    return state

def is_goal(state):
    goals = 0
    for row in state:
        for cell in row:
            if cell in ['$', '%']:
                goals += 1
    return goals == 3

def get_player_positions(state):
    players = []
    for i in range(len(state)):
        for j in range(len(state[i])):
            if state[i][j] in ['*', '%']:
                players.append((i, j))
    return players

def get_valid_moves(state, player_pos):
    moves = []
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    for direction, (dy, dx) in directions.items():
        new_y, new_x = player_pos[0] + dy, player_pos[1] + dx
        
        if state[new_y][new_x] in ['-', 'X']:  # Empty space or goal
            moves.append((direction, player_pos))
        elif state[new_y][new_x] in ['@', '$']:  # Box
            next_y, next_x = new_y + dy, new_x + dx
            if state[next_y][next_x] in ['-', 'X']:  # Can push
                moves.append((direction, player_pos))
    
    return moves

def apply_move(state, move, player_pos):
    direction, pos = move
    new_state = copy.deepcopy(state)
    
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    dy, dx = directions[direction]
    
    y, x = pos
    new_y, new_x = y + dy, x + dx
    
    # Handle player movement
    if new_state[y][x] == '*':
        new_state[y][x] = '-'
    else:  # player on goal
        new_state[y][x] = 'X'
        
    if new_state[new_y][new_x] in ['@', '$']:  # Push box
        next_y, next_x = new_y + dy, new_x + dx
        if new_state[new_y][new_x] == '@':
            if new_state[next_y][next_x] == 'X':
                new_state[next_y][next_x] = '$'
            else:
                new_state[next_y][next_x] = '@'
        else:  # box on goal
            if new_state[next_y][next_x] == 'X':
                new_state[next_y][next_x] = '$'
            else:
                new_state[next_y][next_x] = '@'
    
    # Place player in new position
    if new_state[new_y][new_x] == 'X':
        new_state[new_y][new_x] = '%'
    else:
        new_state[new_y][new_x] = '*'
    
    return new_state

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

def solve_sokoban():
    initial_state = get_initial_state()
    queue = deque([(initial_state, "", None)])
    visited = set()
    
    while queue:
        current_state, path, last_player = queue.popleft()
        state_str = state_to_string(current_state)
        
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        if is_goal(current_state):
            return path
            
        players = get_player_positions(current_state)
        for player_pos in players:
            if last_player is not None and player_pos == last_player:
                continue
                
            valid_moves = get_valid_moves(current_state, player_pos)
            for move in valid_moves:
                new_state = apply_move(current_state, move, player_pos)
                new_path = path + move[0]
                queue.append((new_state, new_path, player_pos))
    
    return None

solution = solve_sokoban()
print(solution)