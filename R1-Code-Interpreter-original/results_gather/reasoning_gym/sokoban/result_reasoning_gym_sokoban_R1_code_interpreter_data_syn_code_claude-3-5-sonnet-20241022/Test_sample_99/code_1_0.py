from collections import deque

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, rows, cols):
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    
    # Find all players
    players = []
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['@', '%']:
                players.append((i, j))
    
    for player_pos in players:
        px, py = player_pos
        for move, dx, dy in directions:
            new_x, new_y = px + dx, py + dy
            
            if not is_valid(new_x, new_y, rows, cols):
                continue
                
            if state[new_x][new_y] == '+':
                continue
                
            new_state = [list(row) for row in state]
            
            # If moving to empty space or goal
            if state[new_x][new_y] in ['-', 'X']:
                # Update player position
                new_state[new_x][new_y] = '%' if state[new_x][new_y] == 'X' else '@'
                new_state[px][py] = 'X' if state[px][py] == '%' else '-'
                next_states.append((move, new_state))
                
            # If moving a box
            elif state[new_x][new_y] in ['$', '@']:
                box_x, box_y = new_x + dx, new_y + dy
                if is_valid(box_x, box_y, rows, cols) and state[box_x][box_y] in ['-', 'X']:
                    # Update box position
                    new_state[box_x][box_y] = '$' if state[box_x][box_y] == 'X' else '@'
                    # Update player position
                    new_state[new_x][new_y] = '%' if state[new_x][new_y] == '$' else '@'
                    new_state[px][py] = 'X' if state[px][py] == '%' else '-'
                    next_states.append((move, new_state))
    
    return next_states

def is_solved(state):
    box_count = 0
    goal_count = 0
    boxes_on_goals = 0
    
    for row in state:
        for cell in row:
            if cell in ['X', '%']:
                goal_count += 1
            if cell == '@':
                box_count += 1
            if cell == '$':
                boxes_on_goals += 1
    
    return boxes_on_goals == 2  # We need exactly 2 boxes on goals

def solve_sokoban(initial_state):
    queue = deque([([], initial_state)])
    seen = set()
    max_moves = 20  # Limit search depth
    
    while queue:
        moves, current_state = queue.popleft()
        
        if len(moves) > max_moves:
            continue
            
        state_tuple = tuple(tuple(row) for row in current_state)
        if state_tuple in seen:
            continue
            
        seen.add(state_tuple)
        
        if is_solved(current_state):
            return ''.join(moves)
            
        for move, next_state in get_next_states(current_state, len(current_state), len(current_state[0])):
            queue.append((moves + [move], next_state))
    
    return None

# Initial state
initial_state = [
    list("++++++++"),
    list("+++$X@%+"),
    list("++$+@@@+"),
    list("+++----+"),
    list("+++$--++"),
    list("+++++X++"),
    list("++++++++"),
    list("++++++++")
]

solution = solve_sokoban(initial_state)
print(solution)