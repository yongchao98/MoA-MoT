from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_moves(state):
    rows = len(state)
    cols = len(state[0])
    
    # Find player position
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if state[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    moves = []
    directions = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    
    for dx, dy, move in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or state[new_x][new_y] == '+':
            continue
            
        if state[new_x][new_y] in ['@', '$']:
            box_new_x, box_new_y = new_x + dx, new_y + dy
            if (not is_valid(box_new_x, box_new_y, rows, cols) or 
                state[box_new_x][box_new_y] in ['+', '@', '$']):
                continue
            
            # Valid push move
            new_state = [list(row) for row in state]
            # Move box
            if state[new_x][new_y] == '@':
                new_state[new_x][new_y] = '-'
            else:  # '$'
                new_state[new_x][new_y] = 'X'
                
            if state[box_new_x][box_new_y] == 'X':
                new_state[box_new_x][box_new_y] = '$'
            else:  # '-'
                new_state[box_new_x][box_new_y] = '@'
            
            # Move player
            if state[player_pos[0]][player_pos[1]] == '*':
                new_state[player_pos[0]][player_pos[1]] = '-'
            else:  # '%'
                new_state[player_pos[0]][player_pos[1]] = 'X'
                
            if state[new_x][new_y] == '@':
                new_state[new_x][new_y] = '*'
            else:  # '$'
                new_state[new_x][new_y] = '%'
                
            moves.append((new_state, move))
        else:
            # Simple move
            new_state = [list(row) for row in state]
            
            if state[player_pos[0]][player_pos[1]] == '*':
                new_state[player_pos[0]][player_pos[1]] = '-'
            else:  # '%'
                new_state[player_pos[0]][player_pos[1]] = 'X'
                
            if state[new_x][new_y] == 'X':
                new_state[new_x][new_y] = '%'
            else:  # '-'
                new_state[new_x][new_y] = '*'
                
            moves.append((new_state, move))
    
    return moves

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    seen = {tuple(map(tuple, initial_state))}
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if solved
        is_solved = True
        for row in current_state:
            for cell in row:
                if cell == '@' or cell == 'X':
                    is_solved = False
                    break
            if not is_solved:
                break
        
        if is_solved:
            return path
        
        for new_state, move in get_moves(current_state):
            state_tuple = tuple(map(tuple, new_state))
            if state_tuple not in seen:
                seen.add(state_tuple)
                queue.append((new_state, path + move))
    
    return None

# Initial state
initial_state = [
    list("++++++"),
    list("+X---+"),
    list("++---+"),
    list("+-@--+"),
    list("+----+"),
    list("+-X@-+"),
    list("+$$--+"),
    list("+*@X-+"),
    list("++++++")
]

solution = solve_sokoban(initial_state)
print(solution)