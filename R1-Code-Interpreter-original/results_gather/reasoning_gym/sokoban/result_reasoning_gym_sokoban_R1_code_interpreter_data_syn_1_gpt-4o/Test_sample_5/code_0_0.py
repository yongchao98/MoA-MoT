from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 7),  # Player position
    ((2, 2), (2, 4), (2, 7)),  # Box positions
)

# Define the goal positions
goal_positions = {(1, 2), (3, 5), (4, 3)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9),
    (1, 0), (1, 9),
    (2, 0), (2, 9),
    (3, 0), (3, 9),
    (4, 0), (4, 5), (4, 9),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9),
}

# Define the directions for movement
directions = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

# Function to check if a position is within bounds and not a wall
def is_valid_position(pos):
    return pos not in walls

# Function to move the player and push a box if necessary
def move(state, direction):
    player_pos, box_positions = state
    dx, dy = directions[direction]
    new_player_pos = (player_pos[0] + dx, player_pos[1] + dy)

    if not is_valid_position(new_player_pos):
        return None

    new_box_positions = set(box_positions)
    if new_player_pos in box_positions:
        new_box_pos = (new_player_pos[0] + dx, new_player_pos[1] + dy)
        if not is_valid_position(new_box_pos) or new_box_pos in box_positions:
            return None
        new_box_positions.remove(new_player_pos)
        new_box_positions.add(new_box_pos)

    return new_player_pos, tuple(sorted(new_box_positions))

# Breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, path = queue.popleft()
        player_pos, box_positions = current_state

        if set(box_positions) == goal_positions:
            return path

        for direction in directions:
            new_state = move(current_state, direction)
            if new_state and new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + direction))

    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)