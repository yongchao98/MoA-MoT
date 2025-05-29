from collections import deque

# Define the initial positions of the knights
initial_positions = {
    'w': [(0, 0), (0, 2)],  # A1, C1
    'B': [(2, 1), (2, 3)]   # B3, D3
}

# Define the target positions for the knights
target_positions = {
    'w': [(2, 1), (2, 3)],  # B3, D3
    'B': [(0, 0), (0, 2)]   # A1, C1
}

# Define the possible moves for a knight
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to coordinates
def pos_to_coord(pos):
    return (ord(pos[0]) - ord('A'), int(pos[1]) - 1)

# Convert coordinates to board positions
def coord_to_pos(coord):
    return chr(coord[0] + ord('A')) + str(coord[1] + 1)

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 4

# Check if the current state matches the target state
def is_target_state(state):
    return state['w'] == target_positions['w'] and state['B'] == target_positions['B']

# Generate the next possible states
def generate_next_states(state, turn):
    next_states = []
    for i, (x, y) in enumerate(state[turn]):
        for dx, dy in knight_moves:
            nx, ny = x + dx, y + dy
            if is_within_board(nx, ny) and (nx, ny) not in state['w'] and (nx, ny) not in state['B']:
                new_state = { 'w': state['w'][:], 'B': state['B'][:] }
                new_state[turn][i] = (nx, ny)
                move = f"{turn},{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"
                next_states.append((new_state, move))
    return next_states

# Perform BFS to find the solution
def bfs_solve():
    initial_state = initial_positions
    queue = deque([(initial_state, [], 'w')])  # (state, moves, turn)
    visited = set()
    visited.add((tuple(initial_state['w']), tuple(initial_state['B'])))
    
    while queue:
        current_state, moves, turn = queue.popleft()
        
        if is_target_state(current_state):
            return moves
        
        next_turn = 'B' if turn == 'w' else 'w'
        for next_state, move in generate_next_states(current_state, turn):
            state_tuple = (tuple(next_state['w']), tuple(next_state['B']))
            if state_tuple not in visited:
                visited.add(state_tuple)
                queue.append((next_state, moves + [move], next_turn))
    
    return "No"

# Find the solution
solution = bfs_solve()
print(solution)