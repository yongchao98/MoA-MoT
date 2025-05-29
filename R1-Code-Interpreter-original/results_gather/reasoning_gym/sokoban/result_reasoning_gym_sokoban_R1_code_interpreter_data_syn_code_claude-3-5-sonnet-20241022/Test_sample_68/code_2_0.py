from collections import deque

def find_player(state):
    for i in range(len(state)):
        for j in range(len(state[i])):
            if state[i][j] in ['*', '%']:
                return (i, j)
    return None

def get_valid_moves(state, pos):
    moves = []
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    rows, cols = len(state), len(state[0])
    x, y = pos

    for direction, (dx, dy) in directions.items():
        new_x, new_y = x + dx, y + dy
        
        # Check if new position is within bounds
        if not (0 <= new_x < rows and 0 <= new_y < cols):
            continue
            
        # If next position is empty or goal
        if state[new_x][new_y] in ['-', 'X']:
            moves.append((direction, (new_x, new_y), None))
            
        # If next position is box
        elif state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if (0 <= box_x < rows and 0 <= box_y < cols and 
                state[box_x][box_y] in ['-', 'X']):
                moves.append((direction, (new_x, new_y), (box_x, box_y)))
                
    return moves

def apply_move(state, player_pos, new_pos, box_move=None):
    new_state = [list(row) for row in state]
    px, py = player_pos
    nx, ny = new_pos
    
    # Handle player's old position
    if state[px][py] == '%':
        new_state[px][py] = 'X'
    else:
        new_state[px][py] = '-'
    
    # Handle player's new position
    if box_move:
        bx, by = box_move
        # Move box
        if state[nx][ny] == '$':
            new_state[nx][ny] = '%'  # Player on goal
            if state[bx][by] == 'X':
                new_state[bx][by] = '$'  # Box on goal
            else:
                new_state[bx][by] = '@'  # Box on floor
        else:
            new_state[nx][ny] = '*'  # Player on floor
            if state[bx][by] == 'X':
                new_state[bx][by] = '$'  # Box on goal
            else:
                new_state[bx][by] = '@'  # Box on floor
    else:
        if state[nx][ny] == 'X':
            new_state[nx][ny] = '%'  # Player on goal
        else:
            new_state[nx][ny] = '*'  # Player on floor
            
    return [''.join(row) for row in new_state]

def is_solved(state):
    goals = 0
    boxes_on_goals = 0
    for row in state:
        for cell in row:
            if cell in ['X', '%']:
                goals += 1
            if cell in ['$', '%']:
                boxes_on_goals += 1
    return goals == boxes_on_goals

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = {tuple(initial_state)}
    max_moves = 100  # Limit search depth
    
    while queue and len(queue[0][1]) < max_moves:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
        
        player_pos = find_player(current_state)
        if not player_pos:
            continue
            
        moves = get_valid_moves(current_state, player_pos)
        for direction, new_pos, box_move in moves:
            new_state = apply_move(current_state, player_pos, new_pos, box_move)
            state_tuple = tuple(new_state)
            
            if state_tuple not in visited:
                visited.add(state_tuple)
                queue.append((new_state, path + direction))
    
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