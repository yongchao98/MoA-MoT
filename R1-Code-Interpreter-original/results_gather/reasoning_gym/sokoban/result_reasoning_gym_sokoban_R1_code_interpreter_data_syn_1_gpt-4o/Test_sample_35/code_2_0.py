from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 1),  # Player position
    frozenset(((1, 2), (2, 1), (2, 2)))  # Box positions
)

# Define the goal positions
goal_positions = frozenset({(3, 1), (3, 2), (4, 3)})

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
    (1, 0), (1, 5),
    (2, 0), (2, 5),
    (3, 0), (3, 3), (3, 5),
    (4, 0), (4, 5),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5)
}

# Define possible moves
moves = {
    'L': (0, -1),
    'R': (0, 1),
    'U': (-1, 0),
    'D': (1, 0)
}

def is_goal_state(box_positions):
    return box_positions == goal_positions

def is_valid_position(position, box_positions):
    return position not in walls and position not in box_positions

def move_player(player_pos, direction, box_positions):
    new_player_pos = (player_pos[0] + direction[0], player_pos[1] + direction[1])
    if not is_valid_position(new_player_pos, box_positions):
        return None, None
    if new_player_pos in box_positions:
        new_box_pos = (new_player_pos[0] + direction[0], new_player_pos[1] + direction[1])
        if not is_valid_position(new_box_pos, box_positions):
            return None, None
        new_box_positions = frozenset(
            new_box_pos if pos == new_player_pos else pos for pos in box_positions
        )
        return new_player_pos, new_box_positions
    return new_player_pos, box_positions

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add(initial_state)

    while queue:
        (player_pos, box_positions), path = queue.popleft()

        if is_goal_state(box_positions):
            return path

        for move, direction in moves.items():
            new_player_pos, new_box_positions = move_player(player_pos, direction, box_positions)
            if new_player_pos is not None:
                new_state = (new_player_pos, new_box_positions)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move))

    return None

solution = solve_sokoban(initial_state)
print(solution)