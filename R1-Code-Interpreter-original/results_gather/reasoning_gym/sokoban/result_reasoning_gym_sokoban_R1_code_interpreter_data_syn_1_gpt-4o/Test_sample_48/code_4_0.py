from collections import deque

# Define the initial state of the puzzle
initial_state = (
    ((1, 4), (2, 5), (3, 3), (5, 4), (6, 1), (7, 3)),  # Box positions
    (7, 2)  # Player position
)

# Define the goal positions
goal_positions = {(2, 2), (3, 1), (7, 4)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
    (1, 0), (1, 1), (1, 2), (1, 3), (1, 5), (1, 6), (1, 7),
    (2, 0), (2, 1), (2, 3), (2, 4), (2, 6), (2, 7),
    (3, 0), (3, 2), (3, 4), (3, 5), (3, 6), (3, 7),
    (4, 0), (4, 1), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7),
    (5, 0), (5, 1), (5, 3), (5, 5), (5, 6), (5, 7),
    (6, 0), (6, 2), (6, 5), (6, 6), (6, 7),
    (7, 0), (7, 1), (7, 5), (7, 6), (7, 7),
    (8, 0), (8, 1), (8, 2), (8, 3), (8, 4), (8, 5), (8, 6), (8, 7)
}

# Define the possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

def is_goal_state(state):
    box_positions, _ = state
    return all(pos in goal_positions for pos in box_positions)

def is_valid_position(position, box_positions):
    return position not in walls and position not in box_positions

def move_player(state, move):
    box_positions, player_position = state
    move_offset = moves[move]
    new_player_position = (player_position[0] + move_offset[0], player_position[1] + move_offset[1])

    if not is_valid_position(new_player_position, box_positions):
        return None

    if new_player_position in box_positions:
        box_index = box_positions.index(new_player_position)
        new_box_position = (new_player_position[0] + move_offset[0], new_player_position[1] + move_offset[1])
        if not is_valid_position(new_box_position, box_positions):
            return None
        new_box_positions = list(box_positions)
        new_box_positions[box_index] = new_box_position
        return (tuple(new_box_positions), new_player_position)
    else:
        return (box_positions, new_player_position)

def solve_sokoban():
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
            new_state = move_player(current_state, move)
            if new_state and new_state not in visited:
                queue.append((new_state, path + move))

    return None

solution = solve_sokoban()
print(solution)