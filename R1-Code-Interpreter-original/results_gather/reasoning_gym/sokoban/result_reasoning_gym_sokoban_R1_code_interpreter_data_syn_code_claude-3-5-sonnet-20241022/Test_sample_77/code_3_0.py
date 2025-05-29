from collections import deque

def get_valid_moves(state, pos):
    moves = []
    rows, cols = len(state), len(state[0])
    x, y = pos
    
    for dx, dy, direction in [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]:
        new_x, new_y = x + dx, y + dy
        
        # Check if move is within bounds
        if 0 <= new_x < rows and 0 <= new_y < cols:
            # Moving to empty space or goal
            if state[new_x][new_y] in ['-', 'X']:
                new_state = [row[:] for row in state]
                # Move player
                if state[x][y] == '*':  # Player on goal
                    new_state[x][y] = 'X'
                else:
                    new_state[x][y] = '-'
                
                if state[new_x][new_y] == 'X':  # Moving to goal
                    new_state[new_x][new_y] = '*'
                else:
                    new_state[new_x][new_y] = '@'
                moves.append((new_state, direction))
            
            # Pushing a box
            elif state[new_x][new_y] == '$':
                push_x, push_y = new_x + dx, new_y + dy
                if (0 <= push_x < rows and 0 <= push_y < cols and 
                    state[push_x][push_y] in ['-', 'X']):
                    new_state = [row[:] for row in state]
                    # Move box
                    new_state[push_x][push_y] = '$'
                    # Move player
                    new_state[new_x][new_y] = '@'
                    if state[x][y] == '*':
                        new_state[x][y] = 'X'
                    else:
                        new_state[x][y] = '-'
                    moves.append((new_state, direction))
    return moves

def find_players(state):
    players = []
    for i in range(len(state)):
        for j in range(len(state[i])):
            if state[i][j] in ['@', '*']:
                players.append((i, j))
    return players

def is_solved(state):
    goals = 0
    boxes_on_goals = 0
    for row in state:
        boxes_on_goals += row.count('$')
        goals += row.count('X')
    return boxes_on_goals == 3 and goals == 0

def solve_sokoban(state):
    initial = [list(row) for row in state]
    queue = deque([(initial, "")])
    seen = {str(initial)}
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
            
        players = find_players(current_state)
        for player in players:
            for new_state, move in get_valid_moves(current_state, player):
                state_str = str(new_state)
                if state_str not in seen:
                    seen.add(state_str)
                    queue.append((new_state, path + move))
    
    return None

# Initial state
initial_state = [
    "++++++++",
    "++$++++",
    "+++X--+",
    "+---$-+",
    "+--@@*-+",
    "+X--$-+",
    "++++++++"
]

solution = solve_sokoban(initial_state)
print(solution)