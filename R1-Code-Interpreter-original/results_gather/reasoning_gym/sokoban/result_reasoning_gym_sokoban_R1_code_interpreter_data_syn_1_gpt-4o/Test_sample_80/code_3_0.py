from collections import deque

# Define the initial state of the puzzle
initial_state = (
    ((1, 4), (2, 3), (3, 2), (3, 3)),  # Box positions
    (3, 4)   # Player position
)

# Define the goal positions
goal_positions = {(1, 5), (2, 2), (3, 6), (4, 1)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8),
    (1, 0), (1, 8),
    (2, 0), (2, 8),
    (3, 0), (3, 8),
    (4, 0), (4, 8),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8)
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

def is_goal_state(box_positions):
    return all(pos in goal_positions for pos in box_positions)

def is_valid_position(position, box_positions):
    return position not in walls and position not in box_positions

def move_player(player_pos, direction):
    return (player_pos[0] + direction[0], player_pos[1] + direction[1])

def move_box(box_pos, direction):
    return (box_pos[0] + direction[0], box_pos[1] + direction[1])

def bfs_solve(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()

    while queue:
        state, path = queue.popleft()
        box_positions, player_pos = state

        if is_goal_state(box_positions):
            return path

        if (box_positions, player_pos) in visited:
            continue

        visited.add((box_positions, player_pos))

        for move, direction in moves.items():
            new_player_pos = move_player(player_pos, direction)

            if not is_valid_position(new_player_pos, box_positions):
                continue

            new_box_positions = list(box_positions)

            if new_player_pos in box_positions:
                box_index = box_positions.index(new_player_pos)
                new_box_pos = move_box(new_player_pos, direction)

                if not is_valid_position(new_box_pos, box_positions):
                    continue

                new_box_positions[box_index] = new_box_pos

            queue.append(((tuple(new_box_positions), new_player_pos), path + move))

    return None

solution = bfs_solve(initial_state)
print(solution)