# Define the initial state of the puzzle
initial_state = [
    "+ + + + + + + + +",
    "+ - X - - - @ % +",
    "+ - - - - - @ @ +",
    "+ - - @ - - + X +",
    "+ + + X $ - + + +",
    "+ + + + + + + + +"
]

# Define the moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Function to find the player position
def find_player(state):
    for i, row in enumerate(state):
        if '%' in row:
            return (i, row.index('%'))
    return None

# Function to check if a move is valid
def is_valid_move(state, player_pos, move):
    if player_pos is None:
        return False
    x, y = player_pos
    dx, dy = moves[move]
    new_x, new_y = x + dx, y + dy
    if state[new_x][new_y] in ['-', 'X']:  # Move to empty or goal
        return True
    if state[new_x][new_y] in ['@', '$']:  # Move to box
        # Check if the box can be pushed
        box_new_x, box_new_y = new_x + dx, new_y + dy
        if state[box_new_x][box_new_y] in ['-', 'X']:  # Box can be pushed to empty or goal
            return True
    return False

# Function to apply a move
def apply_move(state, player_pos, move):
    x, y = player_pos
    dx, dy = moves[move]
    new_x, new_y = x + dx, y + dy
    new_state = [list(row) for row in state]
    
    if new_state[new_x][new_y] in ['-', 'X']:  # Move to empty or goal
        new_state[x][y] = '-' if state[x][y] == '%' else 'X'
        new_state[new_x][new_y] = '%' if state[new_x][new_y] == 'X' else '*'
    elif new_state[new_x][new_y] in ['@', '$']:  # Move to box
        box_new_x, box_new_y = new_x + dx, new_y + dy
        new_state[x][y] = '-' if state[x][y] == '%' else 'X'
        new_state[new_x][new_y] = '%' if state[new_x][new_y] == 'X' else '*'
        new_state[box_new_x][box_new_y] = '$' if state[box_new_x][box_new_y] == 'X' else '@'
    
    return [''.join(row) for row in new_state]

# Function to check if the goal state is reached
def is_goal_state(state):
    for row in state:
        if '@' in row:
            return False
    return True

# Function to solve the puzzle using BFS
from collections import deque

def solve_sokoban(initial_state):
    queue = deque([(initial_state, find_player(initial_state), "")])
    visited = set()
    visited.add(tuple(initial_state))
    
    while queue:
        current_state, player_pos, path = queue.popleft()
        
        if is_goal_state(current_state):
            return path
        
        for move in moves:
            if is_valid_move(current_state, player_pos, move):
                new_state = apply_move(current_state, player_pos, move)
                new_state_tuple = tuple(new_state)
                if new_state_tuple not in visited:
                    visited.add(new_state_tuple)
                    queue.append((new_state, find_player(new_state), path + move))
    
    return None

# Solve the puzzle
solution = solve_sokoban(initial_state)
print(solution)