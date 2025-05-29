from collections import deque

def find_player(state):
    for i in range(len(state)):
        for j in range(len(state[i])):
            if state[i][j] in ['*']:
                return (i, j)
    return None

def get_valid_moves(state, pos):
    moves = []
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    rows, cols = len(state), len(state[0])
    x, y = pos

    for move, (dx, dy) in directions.items():
        new_x, new_y = x + dx, y + dy
        
        if not (0 <= new_x < rows and 0 <= new_y < cols):
            continue
            
        # Simple move to empty space or goal
        if state[new_x][new_y] in ['-', 'X']:
            moves.append((move, (new_x, new_y), None))
            
        # Push box
        elif state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if (0 <= box_x < rows and 0 <= box_y < cols and 
                state[box_x][box_y] in ['-', 'X']):
                moves.append((move, (new_x, new_y), (box_x, box_y)))
                
    return moves

def make_move(state, old_pos, new_pos, box_pos=None):
    new_state = [list(row) for row in state]
    ox, oy = old_pos
    nx, ny = new_pos
    
    # Clear old player position
    new_state[ox][oy] = '-'
    
    # Move box if pushing
    if box_pos:
        bx, by = box_pos
        # Move box
        if new_state[bx][by] == 'X':
            new_state[bx][by] = '$'
        else:
            new_state[bx][by] = '@'
        # Place player
        new_state[nx][ny] = '*'
    else:
        # Just move player
        new_state[nx][ny] = '*'
    
    return [''.join(row) for row in new_state]

def is_goal_state(state):
    box_count = 0
    goal_count = 0
    boxes_on_goals = 0
    
    for row in state:
        for cell in row:
            if cell in ['@', '$']:
                box_count += 1
            if cell in ['X', '$']:
                goal_count += 1
            if cell == '$':
                boxes_on_goals += 1
    
    return box_count == goal_count == boxes_on_goals

def solve_sokoban(initial_state):
    start_pos = find_player(initial_state)
    if not start_pos:
        return None
        
    queue = deque([(initial_state, "", start_pos)])
    visited = {tuple(row for row in initial_state)}
    
    while queue:
        current_state, path, player_pos = queue.popleft()
        
        if is_goal_state(current_state):
            return path
        
        moves = get_valid_moves(current_state, player_pos)
        for direction, new_player_pos, box_pos in moves:
            next_state = make_move(current_state, player_pos, new_player_pos, box_pos)
            state_tuple = tuple(next_state)
            
            if state_tuple not in visited:
                visited.add(state_tuple)
                queue.append((next_state, path + direction, new_player_pos))
    
    return None

# Initial state
initial_state = [
    "+++++++++",
    "+++++X+++",
    "+----@X++",
    "++-@@X$++",
    "++-+--+++",
    "++-$+@$X+",
    "++++-+@+",
    "+++X--@*+",
    "+++++++++"
]

solution = solve_sokoban(initial_state)
print(solution if solution else "No solution found")