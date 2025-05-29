from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_moves(state):
    rows = len(state)
    cols = len(state[0])
    player_pos = None
    
    # Find player position
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
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        if state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if (is_valid(box_x, box_y, rows, cols) and 
                state[box_x][box_y] not in ['+', '@', '$']):
                moves.append((move, (new_x, new_y), (box_x, box_y)))
        elif state[new_x][new_y] in ['-', 'X']:
            moves.append((move, (new_x, new_y), None))
            
    return moves

def apply_move(state, player_pos, box_moves, new_player_pos):
    new_state = [list(row) for row in state]
    px, py = player_pos
    
    # Update player position
    if new_state[px][py] == '%':
        new_state[px][py] = 'X'
    else:
        new_state[px][py] = '-'
        
    npx, npy = new_player_pos
    if new_state[npx][npy] == 'X':
        new_state[npx][npy] = '%'
    else:
        new_state[npx][npy] = '*'
        
    # Update box position if there is one
    if box_moves:
        bx, by = box_moves
        if new_state[npx][npy] == '$':
            new_state[npx][npy] = '@'
        if new_state[bx][by] == 'X':
            new_state[bx][by] = '$'
        else:
            new_state[bx][by] = '@'
            
    return [''.join(row) for row in new_state]

def is_solved(state):
    for row in state:
        for cell in row:
            if cell == '@' or cell == 'X':
                return False
    return True

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "", None)])
    seen = {tuple(initial_state)}
    
    while queue:
        current_state, path, player_pos = queue.popleft()
        
        if is_solved(current_state):
            return path
            
        moves = get_moves(current_state)
        for move_dir, new_pos, box_pos in moves:
            if player_pos is None:
                for i, row in enumerate(current_state):
                    for j, cell in enumerate(row):
                        if cell in ['*', '%']:
                            player_pos = (i, j)
                            break
                    if player_pos:
                        break
                            
            new_state = apply_move(current_state, player_pos, box_pos, new_pos)
            state_tuple = tuple(new_state)
            
            if state_tuple not in seen:
                seen.add(state_tuple)
                queue.append((new_state, path + move_dir, new_pos))
    
    return None

# Initialize the puzzle
initial_state = [
    "++++++",
    "+----+",
    "+@-@-+",
    "+X-@-+",
    "+--X*+",
    "+-X@-+",
    "+$--X+",
    "++++++",
]

solution = solve_sokoban(initial_state)
print(solution if solution else "No solution found")