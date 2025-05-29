from collections import deque
import copy

def is_valid_pos(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_moves(state):
    rows = len(state)
    cols = len(state[0])
    
    # Find player position
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    moves = []
    directions = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    
    for dx, dy, move in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid_pos(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        if state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if (is_valid_pos(box_x, box_y, rows, cols) and 
                state[box_x][box_y] not in ['+', '@', '$']):
                moves.append((move, (new_x, new_y), (box_x, box_y)))
        elif state[new_x][new_y] in ['-', 'X']:
            moves.append((move, (new_x, new_y), None))
            
    return moves

def make_move(state, player_pos, box_move, box_pos):
    new_state = [list(row) for row in state]
    px, py = player_pos
    
    # Update player position
    if new_state[px][py] == '*':
        new_state[px][py] = 'X'
    else:
        new_state[px][py] = '-'
        
    old_px, old_py = None, None
    for i in range(len(state)):
        for j in range(len(state[0])):
            if state[i][j] in ['*', '%']:
                old_px, old_py = i, j
                break
        if old_px is not None:
            break
    
    if state[old_px][old_py] == '%':
        new_state[old_px][old_py] = 'X'
    else:
        new_state[old_px][old_py] = '-'
    
    # Update box position if there is one
    if box_move:
        bx, by = box_pos
        if new_state[bx][by] == 'X':
            new_state[bx][by] = '$'
        else:
            new_state[bx][by] = '@'
            
        if state[px][py] == 'X':
            new_state[px][py] = '%'
        else:
            new_state[px][py] = '*'
    else:
        if state[px][py] == 'X':
            new_state[px][py] = '%'
        else:
            new_state[px][py] = '*'
            
    return [''.join(row) for row in new_state]

def is_solved(state):
    for row in state:
        for cell in row:
            if cell in ['@', 'X']:  # If there's either an unpaired box or goal
                return False
    return True

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = {tuple(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
            
        moves = get_moves(current_state)
        for move, player_pos, box_pos in moves:
            new_state = make_move(current_state, player_pos, box_pos is not None, box_pos)
            state_tuple = tuple(new_state)
            
            if state_tuple not in visited:
                visited.add(state_tuple)
                queue.append((new_state, path + move))
    
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