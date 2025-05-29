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
            new_state = apply_move(state, (x, y), (new_x, new_y))
            moves.append((direction, new_state))
            
        # If next position is box
        elif state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if (0 <= box_x < rows and 0 <= box_y < cols and 
                state[box_x][box_y] in ['-', 'X']):
                new_state = apply_move(state, (x, y), (new_x, new_y), (box_x, box_y))
                moves.append((direction, new_state))
                
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
    
    # Handle player's new position and box movement
    if box_move:
        bx, by = box_move
        # Move box
        if state[nx][ny] == '$':
            new_state[nx][ny] = '%'
        else:
            new_state[nx][ny] = '*'
            
        if state[bx][by] == 'X':
            new_state[bx][by] = '$'
        else:
            new_state[bx][by] = '@'
    else:
        if state[nx][ny] == 'X':
            new_state[nx][ny] = '%'
        else:
            new_state[nx][ny] = '*'
            
    return [''.join(row) for row in new_state]

def is_solved(state):
    for row in state:
        for cell in row:
            if cell == 'X' or cell == '@':  # If there's an unpaired goal or box
                return False
    return True

def is_deadlock(state):
    # Simple corner deadlock detection
    rows, cols = len(state), len(state[0])
    for i in range(rows):
        for j in range(cols):
            if state[i][j] == '@':
                # Check if box is in corner
                if (i > 0 and i < rows-1 and j > 0 and j < cols-1):
                    if ((state[i-1][j] == '+' and state[i][j-1] == '+') or
                        (state[i-1][j] == '+' and state[i][j+1] == '+') or
                        (state[i+1][j] == '+' and state[i][j-1] == '+') or
                        (state[i+1][j] == '+' and state[i][j+1] == '+')):
                        return True
    return False

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = {tuple(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
            
        if is_deadlock(current_state):
            continue
            
        player_pos = find_player(current_state)
        if not player_pos:
            continue
            
        moves = get_valid_moves(current_state, player_pos)
        for direction, new_state in moves:
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