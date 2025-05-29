from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 1),  # Player position
    ((1, 2), (2, 1), (2, 3)),  # Box positions
)

# Define the goal positions
goal_positions = {(1, 2), (2, 4), (4, 4)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
    (1, 0), (1, 6),
    (2, 0), (2, 6),
    (3, 0), (3, 6),
    (4, 0), (4, 2), (4, 6),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6),
}

# Define possible moves
moves = {
    'L': (0, -1),
    'R': (0, 1),
    'U': (-1, 0),
    'D': (1, 0),
}

def is_goal_state(box_positions):
    return all(pos in goal_positions for pos in box_positions)

def is_valid_position(pos, box_positions):
    return pos not in walls and pos not in box_positions

def move_player(state, direction):
    player_pos, box_positions = state
    dx, dy = moves[direction]
    new_player_pos = (player_pos[0] + dx, player_pos[1] + dy)

    if not is_valid_position(new_player_pos, box_positions):
        return None

    if new_player_pos in box_positions:
        new_box_pos = (new_player_pos[0] + dx, new_player_pos[1] + dy)
        if not is_valid_position(new_box_pos, box_positions):
            return None
        new_box_positions = tuple(
            new_box_pos if pos == new_player_pos else pos for pos in box_positions
        )
    else:
        new_box_positions = box_positions

    return new_player_pos, new_box_positions

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()

    while queue:
        current_state, path = queue.popleft()

        if is_goal_state(current_state[1]):
            return path

        if current_state in visited:
            continue

        visited.add(current_state)

        for direction in moves:
            new_state = move_player(current_state, direction)
            if new_state and new_state not in visited:
                queue.append((new_state, path + direction))

    return None

solution = solve_sokoban(initial_state)
print(solution)