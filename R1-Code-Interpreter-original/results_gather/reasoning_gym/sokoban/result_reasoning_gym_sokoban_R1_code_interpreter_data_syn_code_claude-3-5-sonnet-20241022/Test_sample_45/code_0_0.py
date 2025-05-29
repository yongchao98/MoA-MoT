from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_moves(state):
    rows = len(state)
    cols = len(state[0])
    
    # Find player position
    player_x, player_y = None, None
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['*', '%']:
                player_x, player_y = i, j
                break
        if player_x is not None:
            break
    
    moves = []
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    
    for direction, dx, dy in directions:
        new_x, new_y = player_x + dx, player_y + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [list(row) for row in state]
        
        # If moving to empty space or goal
        if state[new_x][new_y] in ['-', 'X']:
            # Update player position
            new_state[new_x][new_y] = '*' if state[new_x][new_y] == '-' else '%'
            new_state[player_x][player_y] = '-' if state[player_x][player_y] == '*' else 'X'
            moves.append((direction, new_state))
            
        # If moving box
        elif state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if (is_valid(box_x, box_y, rows, cols) and 
                state[box_x][box_y] in ['-', 'X']):
                # Update box position
                new_state[box_x][box_y] = '$' if state[box_x][box_y] == 'X' else '@'
                new_state[new_x][new_y] = '%' if state[new_x][new_y] == '$' else '*'
                new_state[player_x][player_y] = '-' if state[player_x][player_y] == '*' else 'X'
                moves.append((direction, new_state))
                
    return moves

def is_solved(state):
    for row in state:
        for cell in row:
            if cell == '@' or cell == 'X':
                return False
    return True

def solve_sokoban(initial_state):
    queue = deque([([], initial_state)])
    visited = set()
    
    while queue:
        path, current_state = queue.popleft()
        state_tuple = tuple(tuple(row) for row in current_state)
        
        if state_tuple in visited:
            continue
            
        visited.add(state_tuple)
        
        if is_solved(current_state):
            return ''.join(path)
            
        for move, new_state in get_moves(current_state):
            queue.append((path + [move], new_state))
            
    return None

# Initial state
initial_state = [
    list("++++++++"),
    list("+++X-+++"),
    list("+X@----+"),
    list("++--X+-+"),
    list("+------+"),
    list("+-+@@$$+"),
    list("+----$*+"),
    list("++++++++")
]

solution = solve_sokoban(initial_state)
print(solution)