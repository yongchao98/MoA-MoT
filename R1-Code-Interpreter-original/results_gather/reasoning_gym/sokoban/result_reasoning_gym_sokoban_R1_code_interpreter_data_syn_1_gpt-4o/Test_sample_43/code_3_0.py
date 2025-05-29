from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 6),  # Player position
    ((2, 4), (2, 5), (3, 4), (4, 2), (4, 5))  # Box positions
)

# Define the goal positions
goal_positions = {(3, 2), (3, 3), (4, 1), (4, 4)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
    (1, 0), (1, 7),
    (2, 0), (2, 7),
    (3, 0), (3, 7),
    (4, 0), (4, 7),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7)
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Function to check if a position is within bounds and not a wall
def is_valid_position(pos):
    return pos not in walls

# Function to apply a move
def apply_move(state, move):
    player_pos, box_positions = state
    move_delta = moves[move]
    new_player_pos = (player_pos[0] + move_delta[0], player_pos[1] + move_delta[1])
    
    if not is_valid_position(new_player_pos):
        return None
    
    new_box_positions = set(box_positions)
    if new_player_pos in box_positions:
        new_box_pos = (new_player_pos[0] + move_delta[0], new_player_pos[1] + move_delta[1])
        if not is_valid_position(new_box_pos) or new_box_pos in box_positions:
            return None
        new_box_positions.remove(new_player_pos)
        new_box_positions.add(new_box_pos)
    
    return (new_player_pos, tuple(sorted(new_box_positions)))

# Breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if set(current_state[1]) == goal_positions:
            return path
        
        for move in moves:
            new_state = apply_move(current_state, move)
            if new_state and new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + move))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)