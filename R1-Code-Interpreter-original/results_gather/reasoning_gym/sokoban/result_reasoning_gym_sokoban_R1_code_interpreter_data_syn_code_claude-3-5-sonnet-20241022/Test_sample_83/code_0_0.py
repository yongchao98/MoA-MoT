from collections import deque
import copy

def is_goal(state):
    return all(cell != '@' for row in state for cell in row)

def get_player_pos(state):
    for i in range(len(state)):
        for j in range(len(state[i])):
            if state[i][j] in ['*', '%']:
                return (i, j)
    return None

def is_valid_move(state, player_pos, dx, dy):
    new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
    if state[new_x][new_y] == '+':
        return False
    if state[new_x][new_y] in ['@', '$']:
        box_x, box_y = new_x + dx, new_y + dy
        if state[box_x][box_y] in ['+', '@', '$']:
            return False
    return True

def make_move(state, player_pos, dx, dy):
    new_state = [list(row) for row in state]
    x, y = player_pos
    new_x, new_y = x + dx, y + dy
    
    # Handle player movement
    if new_state[x][y] == '*':
        new_state[x][y] = '-'
    elif new_state[x][y] == '%':
        new_state[x][y] = 'X'
        
    if new_state[new_x][new_y] == '-':
        new_state[new_x][new_y] = '*'
    elif new_state[new_x][new_y] == 'X':
        new_state[new_x][new_y] = '%'
    elif new_state[new_x][new_y] in ['@', '$']:
        box_x, box_y = new_x + dx, new_y + dy
        if new_state[new_x][new_y] == '@':
            if new_state[box_x][box_y] == 'X':
                new_state[box_x][box_y] = '$'
            else:
                new_state[box_x][box_y] = '@'
        else:  # '$'
            if new_state[box_x][box_y] == 'X':
                new_state[box_x][box_y] = '$'
            else:
                new_state[box_x][box_y] = '@'
        if new_state[new_x][new_y] == '@':
            new_state[new_x][new_y] = '*'
        else:  # '$'
            new_state[new_x][new_y] = '%'
            
    return [''.join(row) for row in new_state]

def solve_sokoban(initial_state):
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    queue = deque([(initial_state, "", get_player_pos(initial_state))])
    visited = {tuple(initial_state)}
    
    while queue:
        current_state, path, player_pos = queue.popleft()
        
        if is_goal(current_state):
            return path
            
        for move, (dx, dy) in moves.items():
            if is_valid_move(current_state, player_pos, dx, dy):
                new_state = make_move(current_state, player_pos, dx, dy)
                new_state_tuple = tuple(new_state)
                
                if new_state_tuple not in visited:
                    visited.add(new_state_tuple)
                    new_player_pos = (player_pos[0] + dx, player_pos[1] + dy)
                    queue.append((new_state, path + move, new_player_pos))
    
    return None

initial_state = [
    "++++++++",
    "+++$*$++",
    "+---@-++",
    "+-@-X-$+",
    "+--+++++",
    "+--X++++",
    "++++++++",
]

solution = solve_sokoban(initial_state)
print(solution if solution else "No solution found")