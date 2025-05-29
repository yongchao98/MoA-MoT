from collections import deque

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, player_pos):
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    rows, cols = len(state), len(state[0])
    
    for move, dx, dy in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [row[:] for row in state]
        new_player_pos = [new_x, new_y]
        
        # Moving to empty space or goal
        if state[new_x][new_y] in '-X':
            old_cell = 'X' if state[player_pos[0]][player_pos[1]] == '%' else '-'
            new_cell = '%' if state[new_x][new_y] == 'X' else '*'
            new_state[player_pos[0]][player_pos[1]] = old_cell
            new_state[new_x][new_y] = new_cell
            next_states.append((move, new_state, new_player_pos))
        
        # Moving a box
        elif state[new_x][new_y] in '@$':
            box_x, box_y = new_x + dx, new_y + dy
            if is_valid(box_x, box_y, rows, cols) and state[box_x][box_y] in '-X':
                old_cell = 'X' if state[player_pos[0]][player_pos[1]] == '%' else '-'
                mid_cell = '%' if state[new_x][new_y] == '$' else '*'
                new_cell = '$' if state[box_x][box_y] == 'X' else '@'
                
                new_state[player_pos[0]][player_pos[1]] = old_cell
                new_state[new_x][new_y] = mid_cell
                new_state[box_x][box_y] = new_cell
                next_states.append((move, new_state, new_player_pos))
    
    return next_states

def state_to_string(state):
    return '\n'.join(''.join(row) for row in state)

def is_solved(state):
    goals = 0
    boxes_on_goals = 0
    for row in state:
        for cell in row:
            if cell in 'X':
                goals += 1
            elif cell in '$':
                boxes_on_goals += 1
    return goals == 0 and boxes_on_goals > 0

def solve_sokoban():
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
    
    # Find player position
    player_pos = None
    for i in range(len(initial_state)):
        for j in range(len(initial_state[i])):
            if initial_state[i][j] in '*%':
                player_pos = [i, j]
                break
        if player_pos:
            break
    
    queue = deque([([], initial_state, player_pos)])
    visited = {state_to_string(initial_state)}
    
    while queue:
        moves, current_state, current_pos = queue.popleft()
        
        if len(moves) > 30:  # Limit search depth
            continue
            
        if is_solved(current_state):
            return ''.join(moves)
        
        for move, next_state, next_pos in get_next_states(current_state, current_pos):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((moves + [move], next_state, next_pos))
    
    return "No solution found"

print(solve_sokoban())