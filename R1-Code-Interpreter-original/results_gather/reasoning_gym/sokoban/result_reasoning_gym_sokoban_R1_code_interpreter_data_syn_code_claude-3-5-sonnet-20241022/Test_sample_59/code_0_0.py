from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, rows, cols):
    directions = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    next_states = []
    
    # Find player position
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    for dx, dy, move in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        new_state = [list(row) for row in state]
        
        # If moving to empty space or goal
        if state[new_x][new_y] in ['-', 'X']:
            new_state[new_x][new_y] = '*' if state[new_x][new_y] == '-' else '%'
            new_state[player_pos[0]][player_pos[1]] = '-' if state[player_pos[0]][player_pos[1]] == '*' else 'X'
            next_states.append((new_state, move))
            
        # If pushing a box
        elif state[new_x][new_y] in ['@', '$']:
            push_x, push_y = new_x + dx, new_y + dy
            if (is_valid(push_x, push_y, rows, cols) and 
                state[push_x][push_y] in ['-', 'X']):
                new_state[push_x][push_y] = '$' if state[push_x][push_y] == 'X' else '@'
                new_state[new_x][new_y] = '*' if state[new_x][new_y] == '@' else '%'
                new_state[player_pos[0]][player_pos[1]] = '-' if state[player_pos[0]][player_pos[1]] == '*' else 'X'
                next_states.append((new_state, move))
                
    return next_states

def is_solved(state):
    return all(box == '$' for row in state for box in row if box in ['@', '$'])

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

def solve_sokoban(initial_state):
    rows, cols = len(initial_state), len(initial_state[0])
    queue = deque([(initial_state, "")])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if is_solved(current_state):
            return moves
            
        for next_state, move in get_next_states(current_state, rows, cols):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, moves + move))
    
    return None

# Initial state
initial_state = [
    list("++++++++"),
    list("+-X@-XX+"),
    list("+----@-+"),
    list("+---X--+"),
    list("+-@@--++"),
    list("+-@---X+"),
    list("+X@-$--+"),
    list("+*-----+"),
    list("++++++++")
]

solution = solve_sokoban(initial_state)
print(solution)