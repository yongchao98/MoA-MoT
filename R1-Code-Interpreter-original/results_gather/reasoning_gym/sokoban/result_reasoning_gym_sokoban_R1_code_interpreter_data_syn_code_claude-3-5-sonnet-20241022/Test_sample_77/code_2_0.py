from collections import deque

def get_state_string(state):
    return '\n'.join(''.join(row) for row in state)

def get_moves(state):
    rows = len(state)
    cols = len(state[0])
    moves = []
    
    # Find player positions
    players = []
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['@', '*']:
                players.append((i, j))

    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    for player_pos in players:
        px, py = player_pos
        for move, (dx, dy) in directions.items():
            new_x, new_y = px + dx, py + dy
            
            # Check if move is valid
            if (0 <= new_x < rows and 0 <= new_y < cols and 
                state[new_x][new_y] not in ['+', '@', '*']):
                
                # Create new state
                new_state = [list(row) for row in state]
                
                # Handle player movement
                if state[new_x][new_y] in ['-', 'X']:
                    # Move player
                    if state[px][py] == '*':  # Player on goal
                        new_state[px][py] = 'X'
                    else:
                        new_state[px][py] = '-'
                    
                    if state[new_x][new_y] == 'X':  # Moving to goal
                        new_state[new_x][new_y] = '*'
                    else:
                        new_state[new_x][new_y] = '@'
                    
                    moves.append((new_state, move))
                
                # Handle box pushing
                elif state[new_x][new_y] == '$':
                    push_x, push_y = new_x + dx, new_y + dy
                    
                    if (0 <= push_x < rows and 0 <= push_y < cols and 
                        state[push_x][push_y] in ['-', 'X']):
                        
                        # Move box
                        if state[push_x][push_y] == 'X':
                            new_state[push_x][push_y] = '$'
                        else:
                            new_state[push_x][push_y] = '$'
                            
                        # Move player
                        new_state[new_x][new_y] = '@'
                        if state[px][py] == '*':
                            new_state[px][py] = 'X'
                        else:
                            new_state[px][py] = '-'
                            
                        moves.append((new_state, move))
    
    return moves

def is_solved(state):
    box_count = 0
    goal_count = 0
    for row in state:
        box_count += row.count('$')
        goal_count += row.count('X')
    return box_count == 3 and goal_count == 0

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    seen = {get_state_string(initial_state)}
    max_moves = 20  # Limit search depth
    
    while queue:
        current_state, path = queue.popleft()
        
        if len(path) > max_moves:
            continue
            
        if is_solved(current_state):
            return path
        
        for new_state, move in get_moves(current_state):
            state_string = get_state_string(new_state)
            if state_string not in seen:
                seen.add(state_string)
                queue.append((new_state, path + move))
    
    return None

# Initial state
initial_state = [
    list("++++++++"),
    list("++$++++"),
    list("+++X--+"),
    list("+---$-+"),
    list("+--@@*-+"),
    list("+X--$-+"),
    list("++++++++")
]

solution = solve_sokoban(initial_state)
print(solution)