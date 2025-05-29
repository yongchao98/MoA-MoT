from collections import deque

# Define the initial state of the puzzle
initial_state = [
    list("+++++++"),
    list("++++-$+"),
    list("++----+"),
    list("+-----+"),
    list("+X----+"),
    list("+---@X+"),
    list("+--@+-+"),
    list("+$$X-@+"),
    list("++++$*+"),
    list("+++++++")
]

# Define the possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Function to find the player's position
def find_player(state):
    for i, row in enumerate(state):
        if '*' in row:
            return i, row.index('*')
    return None

# Function to check if a move is valid
def is_valid_move(state, move, player_pos):
    x, y = player_pos
    dx, dy = moves[move]
    new_x, new_y = x + dx, y + dy
    # Check boundaries
    if new_x < 0 or new_x >= len(state) or new_y < 0 or new_y >= len(state[0]):
        return False
    if state[new_x][new_y] in ('+', '@', '$'):
        return False
    return True

# Function to apply a move
def apply_move(state, move, player_pos):
    x, y = player_pos
    dx, dy = moves[move]
    new_x, new_y = x + dx, y + dy
    new_state = [row[:] for row in state]
    new_state[x][y] = '-'
    new_state[new_x][new_y] = '*'
    return new_state, (new_x, new_y)

# Breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, find_player(initial_state), "")])
    visited = set()
    visited.add(tuple(map(tuple, initial_state)))
    
    while queue:
        current_state, player_pos, path = queue.popleft()
        
        # Check if all boxes are on goals
        if all(row.count('@') == 0 for row in current_state):
            return path
        
        for move in moves:
            if is_valid_move(current_state, move, player_pos):
                new_state, new_player_pos = apply_move(current_state, move, player_pos)
                state_tuple = tuple(map(tuple, new_state))
                if state_tuple not in visited:
                    visited.add(state_tuple)
                    queue.append((new_state, new_player_pos, path + move))
    
    return "No solution"

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)