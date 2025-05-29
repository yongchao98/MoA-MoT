from collections import deque
import copy

def get_initial_state():
    puzzle = [
        "+ + + + + + + + + +",
        "+ - - X - - - @ * +",
        "+ - - - @ X @ @ @ +",
        "+ - - - X + X $ - +",
        "+ + - - + $ + + X +",
        "+ X - - + + + + $ +",
        "+ - @ - + + + + + +",
        "+ - - - - + + + + +",
        "+ - - - - + + + $ +",
        "+ + + + + + + + + +"
    ]
    return [list(row.replace(" ", "")) for row in puzzle]

def count_goals_and_boxes(state):
    goals = 0
    boxes = 0
    boxes_on_goals = 0
    for row in state:
        for cell in row:
            if cell in ['X', '*', '$']:  # All goal positions
                goals += 1
            if cell in ['@', '*', '$']:  # All boxes
                boxes += 1
            if cell in ['$', '*']:  # Boxes on goals
                boxes_on_goals += 1
    return goals, boxes, boxes_on_goals

def is_goal_state(state):
    goals, boxes, boxes_on_goals = count_goals_and_boxes(state)
    return boxes == boxes_on_goals == goals

def get_player_pos(state):
    for i in range(len(state)):
        for j in range(len(state[i])):
            if state[i][j] in ['@', '*']:
                return (i, j)
    return None

def is_deadlock(state, r, c):
    # Simple corner deadlock detection
    if state[r][c] in ['@', '$']:
        corners = [(r-1, c), (r+1, c), (r, c-1), (r, c+1)]
        wall_count = sum(1 for nr, nc in corners if 0 <= nr < len(state) and 0 <= nc < len(state[0]) and state[nr][nc] == '+')
        if wall_count >= 2:
            return True
    return False

def get_valid_moves(state):
    moves = []
    player_pos = get_player_pos(state)
    if not player_pos:
        return moves
    
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    for direction, (dr, dc) in directions.items():
        new_state = [row[:] for row in state]
        r, c = player_pos
        new_r, new_c = r + dr, c + dc
        
        if 0 <= new_r < len(state) and 0 <= new_c < len(state[0]):
            curr_char = state[r][c]
            next_char = state[new_r][new_c]
            
            # Moving to empty space or goal
            if next_char in ['-', 'X']:
                new_state[r][c] = 'X' if curr_char == '*' else '-'
                new_state[new_r][new_c] = '*' if next_char == 'X' else '@'
                moves.append((direction, new_state))
            
            # Pushing a box
            elif next_char in ['@', '$']:
                push_r, push_c = new_r + dr, new_c + dc
                if 0 <= push_r < len(state) and 0 <= push_c < len(state[0]):
                    push_char = state[push_r][push_c]
                    if push_char in ['-', 'X']:
                        # Check for deadlock before pushing
                        if not is_deadlock(state, push_r, push_c):
                            new_state[r][c] = 'X' if curr_char == '*' else '-'
                            new_state[new_r][new_c] = '*' if next_char == '$' else '@'
                            new_state[push_r][push_c] = '$' if push_char == 'X' else '@'
                            moves.append((direction, new_state))
    
    return moves

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

def solve_sokoban():
    initial_state = get_initial_state()
    queue = deque([(initial_state, "")])
    visited = {state_to_string(initial_state)}
    max_moves = 100  # Limit the search depth
    
    while queue:
        state, path = queue.popleft()
        
        if len(path) > max_moves:
            continue
            
        if is_goal_state(state):
            return path
        
        for move, new_state in get_valid_moves(state):
            state_str = state_to_string(new_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((new_state, path + move))
    
    return None

solution = solve_sokoban()
print(solution if solution else "No solution found")