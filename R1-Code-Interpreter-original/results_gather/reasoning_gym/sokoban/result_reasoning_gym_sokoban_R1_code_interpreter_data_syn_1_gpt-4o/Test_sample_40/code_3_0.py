from collections import deque

# Define the initial state
initial_state = (
    (1, 1),  # Player position
    frozenset({(1, 2), (2, 1), (2, 2), (2, 3), (3, 3)}),  # Box positions
)

# Define the goal positions
goal_positions = frozenset({(1, 4), (3, 1), (6, 3)})

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
    (1, 0), (1, 5),
    (2, 0), (2, 5), (2, 6),
    (3, 0), (3, 5), (3, 6),
    (4, 0), (4, 3), (4, 5),
    (5, 0), (5, 3), (5, 4), (5, 5),
    (6, 0), (6, 3), (6, 4), (6, 5),
    (7, 0), (7, 2), (7, 3), (7, 4), (7, 5),
    (8, 0), (8, 1), (8, 2), (8, 3), (8, 4), (8, 5),
    (9, 0), (9, 1), (9, 2), (9, 3), (9, 4), (9, 5),
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

def is_goal_state(state):
    _, box_positions = state
    return box_positions == goal_positions

def is_valid_position(pos, box_positions):
    return pos not in walls and pos not in box_positions

def move(state, direction):
    player_pos, box_positions = state
    dx, dy = moves[direction]
    new_player_pos = (player_pos[0] + dx, player_pos[1] + dy)

    if not is_valid_position(new_player_pos, box_positions):
        return None

    new_box_positions = set(box_positions)
    if new_player_pos in box_positions:
        new_box_pos = (new_player_pos[0] + dx, new_player_pos[1] + dy)
        if not is_valid_position(new_box_pos, box_positions):
            return None
        new_box_positions.remove(new_player_pos)
        new_box_positions.add(new_box_pos)

    return (new_player_pos, frozenset(new_box_positions))

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

        for direction in moves:
            new_state = move(current_state, direction)
            if new_state and new_state not in visited:
                queue.append((new_state, path + direction))

    return None

solution = solve_sokoban(initial_state)
print(solution)