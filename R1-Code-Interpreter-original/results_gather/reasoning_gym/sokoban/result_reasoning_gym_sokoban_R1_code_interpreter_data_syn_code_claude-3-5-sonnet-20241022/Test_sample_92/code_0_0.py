from collections import deque
import copy

# Initial state representation
initial_state = [
    "++++++++",
    "+++X@-$+",
    "+*--@--+",
    "+-@X--X+",
    "+-$----+",
    "++++++++"
]

def is_goal(state):
    # Check if all boxes are on goals
    boxes = 0
    goals = 0
    boxes_on_goals = 0
    for row in state:
        for cell in row:
            if cell in '@$':
                boxes += 1
            if cell in 'X$':
                goals += 1
            if cell == '$':
                boxes_on_goals += 1
    return boxes == goals == boxes_on_goals

def get_player_pos(state):
    for i in range(len(state)):
        for j in range(len(state[i])):
            if state[i][j] in '@%':
                return (i, j)
    return None

def is_valid_move(state, player_pos, dx, dy):
    new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
    if state[new_x][new_y] == '+':
        return False
    if state[new_x][new_y] in '@$':
        box_x, box_y = new_x + dx, new_y + dy
        if state[box_x][box_y] in '+@$':
            return False
    return True

def make_move(state, direction):
    state = [list(row) for row in state]
    player_pos = get_player_pos(state)
    if not player_pos:
        return None
    
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    dx, dy = moves[direction]
    
    if not is_valid_move(state, player_pos, dx, dy):
        return None
        
    px, py = player_pos
    new_px, new_py = px + dx, py + dy
    
    # Handle box pushing
    if state[new_px][new_py] in '@$':
        box_x, box_y = new_px + dx, new_py + dy
        # Move box
        state[box_x][box_y] = '$' if state[box_x][box_y] == 'X' else '@'
        
    # Update player position
    state[new_px][new_py] = '%' if state[new_px][new_py] in 'X$' else '@'
    state[px][py] = 'X' if state[px][py] == '%' else '-'
    
    return [''.join(row) for row in state]

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set([''.join(initial_state)])
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal(current_state):
            return path
            
        for move in 'UDLR':
            new_state = make_move(current_state, move)
            if new_state:
                state_str = ''.join(new_state)
                if state_str not in visited:
                    visited.add(state_str)
                    queue.append((new_state, path + move))
    
    return None

solution = solve_sokoban(initial_state)
print(solution)