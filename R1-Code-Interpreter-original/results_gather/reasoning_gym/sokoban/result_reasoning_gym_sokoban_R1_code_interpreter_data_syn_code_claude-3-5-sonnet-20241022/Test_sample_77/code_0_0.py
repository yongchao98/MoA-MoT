from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_moves(state):
    rows = len(state)
    cols = len(state[0])
    moves = []
    
    # Find players (@ and *)
    players = []
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['@', '*']:
                players.append((i, j))
    
    directions = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    
    for player in players:
        px, py = player
        for dx, dy, move in directions:
            new_x, new_y = px + dx, py + dy
            
            if not is_valid(new_x, new_y, rows, cols):
                continue
                
            if state[new_x][new_y] == '+':
                continue
                
            # If moving to empty space or goal
            if state[new_x][new_y] in ['-', 'X']:
                new_state = [list(row) for row in state]
                # Update player position
                if state[px][py] == '@':
                    new_state[px][py] = '-'
                else:  # player is on goal (*)
                    new_state[px][py] = 'X'
                    
                if state[new_x][new_y] == '-':
                    new_state[new_x][new_y] = '@'
                else:  # moving to goal
                    new_state[new_x][new_y] = '*'
                    
                moves.append(([''.join(row) for row in new_state], move))
                
            # If pushing a box
            elif state[new_x][new_y] in ['$', '@']:
                push_x, push_y = new_x + dx, new_y + dy
                
                if not is_valid(push_x, push_y, rows, cols):
                    continue
                    
                if state[push_x][push_y] in ['-', 'X']:
                    new_state = [list(row) for row in state]
                    # Update player position
                    if state[px][py] == '@':
                        new_state[px][py] = '-'
                    else:  # player is on goal (*)
                        new_state[px][py] = 'X'
                        
                    # Update box position
                    if state[push_x][push_y] == '-':
                        new_state[push_x][push_y] = '$'
                    else:  # box pushed to goal
                        new_state[push_x][push_y] = '$'
                        
                    # Update intermediate position
                    if state[new_x][new_y] == '$':
                        new_state[new_x][new_y] = '@'
                    else:
                        new_state[new_x][new_y] = '@'
                        
                    moves.append(([''.join(row) for row in new_state], move))
    
    return moves

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    seen = {tuple(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if solved
        goals_covered = 0
        for row in current_state:
            goals_covered += row.count('$')
        if goals_covered == 3:  # number of boxes needed on goals
            return path
        
        # Try all possible moves
        for new_state, move in get_moves(current_state):
            state_tuple = tuple(new_state)
            if state_tuple not in seen:
                seen.add(state_tuple)
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
    "++++++++",
]

solution = solve_sokoban(initial_state)
print(solution)