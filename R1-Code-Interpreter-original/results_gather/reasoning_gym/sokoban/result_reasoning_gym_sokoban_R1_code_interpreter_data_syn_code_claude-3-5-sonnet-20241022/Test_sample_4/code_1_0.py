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
        
        # Handle box pushing
        if state[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if not is_valid(box_x, box_y, rows, cols) or state[box_x][box_y] in ['+', '@', '$']:
                continue
            
            # Move box
            new_state[box_x][box_y] = '$' if state[box_x][box_y] == 'X' else '@'
            new_state[new_x][new_y] = 'X' if state[new_x][new_y] == '$' else '-'
        
        # Move player
        new_state[player_pos[0]][player_pos[1]] = 'X' if state[player_pos[0]][player_pos[1]] == '*' else '-'
        new_state[new_x][new_y] = '*' if state[new_x][new_y] == 'X' else '*'
        
        next_states.append((new_state, (new_x, new_y), direction))
    return next_states

def state_to_string(state):
    return '\n'.join(''.join(row) for row in state)

def count_goals_and_boxes_on_goals(state):
    goals = 0
    boxes_on_goals = 0
    for row in state:
        for cell in row:
            if cell in ['X', '*']:
                goals += 1
            if cell == '$':
                boxes_on_goals += 1
    return goals, boxes_on_goals

def solve_sokoban(initial_state):
    rows, cols = len(initial_state), len(initial_state[0])
    total_goals, _ = count_goals_and_boxes_on_goals(initial_state)
    
    # Find player position
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if initial_state[i][j] == '*':
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    queue = deque([(initial_state, player_pos, "")])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, current_pos, path = queue.popleft()
        
        # Check if solved
        _, boxes_on_goals = count_goals_and_boxes_on_goals(current_state)
        if boxes_on_goals == total_goals:
            return path
        
        for next_state, next_pos, direction in get_next_states(current_state, current_pos):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, next_pos, path + direction))
    
    return None

# Initial state
initial_state = [
    list("++++++++"),
    list("+X+----+"),
    list("+*$--@@+"),
    list("+-@@+XX+"),
    list("+--X++++"),
    list("++++++++")
]

solution = solve_sokoban(initial_state)
print(solution)