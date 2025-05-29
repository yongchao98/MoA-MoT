from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, player_pos):
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    rows, cols = len(state), len(state[0])
    
    for direction, dx, dy in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [list(row) for row in state]
        new_player_pos = [new_x, new_y]
        
        # If moving to empty space or goal
        if state[new_x][new_y] in '-X':
            new_state[player_pos[0]][player_pos[1]] = '-' if state[player_pos[0]][player_pos[1]] == '*' else 'X'
            new_state[new_x][new_y] = '*' if state[new_x][new_y] == '-' else '%'
            next_states.append((direction, new_state, new_player_pos))
            
        # If moving a box
        elif state[new_x][new_y] in '@$':
            box_x, box_y = new_x + dx, new_y + dy
            if is_valid(box_x, box_y, rows, cols) and state[box_x][box_y] in '-X':
                new_state[player_pos[0]][player_pos[1]] = '-' if state[player_pos[0]][player_pos[1]] == '*' else 'X'
                new_state[new_x][new_y] = '*' if state[new_x][new_y] == '@' else '%'
                new_state[box_x][box_y] = '@' if state[box_x][box_y] == '-' else '$'
                next_states.append((direction, new_state, new_player_pos))
                
    return next_states

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

def is_solved(state):
    for row in state:
        for cell in row:
            if cell == '@' or cell == 'X':
                return False
    return True

def solve_sokoban(initial_state):
    rows, cols = len(initial_state), len(initial_state[0])
    # Find player position
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if initial_state[i][j] in '*%':
                player_pos = [i, j]
                break
        if player_pos:
            break
    
    queue = deque([([], initial_state, player_pos)])
    visited = {state_to_string(initial_state)}
    
    while queue:
        moves, current_state, current_pos = queue.popleft()
        
        if is_solved(current_state):
            return ''.join(moves)
            
        for direction, next_state, next_pos in get_next_states(current_state, current_pos):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((moves + [direction], next_state, next_pos))
    
    return None

# Initial state
initial_state = [
    list("++++++++"),
    list("+X-X*-+"),
    list("+-@@@@-+"),
    list("+-X--+"),
    list("+---$-+"),
    list("+-X--+"),
    list("++---+"),
    list("+$+X@-+"),
    list("++++++++")
]

solution = solve_sokoban(initial_state)
print(solution)