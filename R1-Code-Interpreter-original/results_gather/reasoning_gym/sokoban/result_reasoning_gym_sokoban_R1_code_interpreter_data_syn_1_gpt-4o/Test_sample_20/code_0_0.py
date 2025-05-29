from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 3), (1, 5), (2, 5),  # Box positions
    (4, 4)                   # Player position
)

# Define the goal positions
goal_positions = {(1, 1), (1, 2), (2, 3), (2, 4)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
    (1, 0), (1, 7),
    (2, 0), (2, 7),
    (3, 0), (3, 7),
    (4, 0), (4, 7),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7)
}

# Define the possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Check if the state is a goal state
def is_goal_state(state):
    box_positions = state[:3]
    return all(pos in goal_positions for pos in box_positions)

# Check if a move is valid
def is_valid_move(state, move):
    box_positions = set(state[:3])
    player_position = state[3]
    move_offset = moves[move]
    new_player_position = (player_position[0] + move_offset[0], player_position[1] + move_offset[1])

    if new_player_position in walls:
        return False

    if new_player_position in box_positions:
        new_box_position = (new_player_position[0] + move_offset[0], new_player_position[1] + move_offset[1])
        if new_box_position in walls or new_box_position in box_positions:
            return False

    return True

# Apply a move to the state
def apply_move(state, move):
    box_positions = list(state[:3])
    player_position = state[3]
    move_offset = moves[move]
    new_player_position = (player_position[0] + move_offset[0], player_position[1] + move_offset[1])

    if new_player_position in box_positions:
        box_index = box_positions.index(new_player_position)
        new_box_position = (new_player_position[0] + move_offset[0], new_player_position[1] + move_offset[1])
        box_positions[box_index] = new_box_position

    return tuple(box_positions + [new_player_position])

# Perform a breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()

    while queue:
        current_state, path = queue.popleft()

        if is_goal_state(current_state):
            return path

        if current_state in visited:
            continue

        visited.add(current_state)

        for move in moves:
            if is_valid_move(current_state, move):
                new_state = apply_move(current_state, move)
                queue.append((new_state, path + move))

    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)