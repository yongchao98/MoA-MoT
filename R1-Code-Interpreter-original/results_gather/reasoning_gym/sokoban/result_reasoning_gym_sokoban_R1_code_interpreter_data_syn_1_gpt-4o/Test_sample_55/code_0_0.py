from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 5),  # Player position
    frozenset([(1, 4), (2, 4), (2, 5), (4, 2), (5, 3)])  # Box positions
)

# Define the goal positions
goal_positions = frozenset([(1, 3), (3, 3), (4, 5)])

# Define the walls
walls = frozenset([
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
    (1, 0), (1, 1), (1, 2), (1, 6),
    (2, 0), (2, 1), (2, 2), (2, 6),
    (3, 0), (3, 1), (3, 2), (3, 6),
    (4, 0), (4, 1), (4, 6),
    (5, 0), (5, 1), (5, 4), (5, 6),
    (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6)
])

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

def is_goal_state(box_positions):
    return box_positions == goal_positions

def move_player(player_pos, direction):
    return (player_pos[0] + direction[0], player_pos[1] + direction[1])

def move_box(box_pos, direction):
    return (box_pos[0] + direction[0], box_pos[1] + direction[1])

def is_valid_position(pos, box_positions):
    return pos not in walls and pos not in box_positions

def solve_sokoban():
    queue = deque([(initial_state, "")])
    visited = set([initial_state])

    while queue:
        (player_pos, box_positions), path = queue.popleft()

        if is_goal_state(box_positions):
            return path

        for move, direction in moves.items():
            new_player_pos = move_player(player_pos, direction)

            if is_valid_position(new_player_pos, box_positions):
                new_box_positions = box_positions

                if new_player_pos in box_positions:
                    new_box_pos = move_box(new_player_pos, direction)

                    if is_valid_position(new_box_pos, box_positions):
                        new_box_positions = frozenset(
                            (new_box_pos if box == new_player_pos else box) for box in box_positions
                        )
                    else:
                        continue

                new_state = (new_player_pos, new_box_positions)

                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move))

    return "No solution found"

solution = solve_sokoban()
print(solution)