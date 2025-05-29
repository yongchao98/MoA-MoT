from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 1),  # Player position
    frozenset([(2, 2), (2, 3), (4, 2)]),  # Box positions
)

# Define the goal positions
goal_positions = frozenset([(2, 1), (3, 1), (4, 3)])

# Define the walls
walls = frozenset([
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
    (1, 0), (1, 6),
    (2, 0), (2, 6),
    (3, 0), (3, 6),
    (4, 0), (4, 6),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6),
])

# Define possible moves
moves = {
    'L': (0, -1),
    'R': (0, 1),
    'U': (-1, 0),
    'D': (1, 0),
}

def is_valid_position(position, boxes):
    """Check if a position is valid (not a wall or a box)."""
    return position not in walls and position not in boxes

def move_player(player, direction):
    """Move the player in the given direction."""
    return (player[0] + direction[0], player[1] + direction[1])

def move_box(box, direction):
    """Move the box in the given direction."""
    return (box[0] + direction[0], box[1] + direction[1])

def solve_sokoban(initial_state, goal_positions):
    """Solve the Sokoban puzzle using BFS."""
    queue = deque([(initial_state, "")])
    visited = set([initial_state])

    while queue:
        (player, boxes), path = queue.popleft()

        if boxes == goal_positions:
            return path

        for move, direction in moves.items():
            new_player = move_player(player, direction)

            if is_valid_position(new_player, boxes):
                if new_player in boxes:
                    new_box = move_box(new_player, direction)
                    if is_valid_position(new_box, boxes):
                        new_boxes = frozenset((new_box if box == new_player else box) for box in boxes)
                        new_state = (new_player, new_boxes)
                        if new_state not in visited:
                            visited.add(new_state)
                            queue.append((new_state, path + move))
                else:
                    new_state = (new_player, boxes)
                    if new_state not in visited:
                        visited.add(new_state)
                        queue.append((new_state, path + move))

    return None

solution = solve_sokoban(initial_state, goal_positions)
print(solution)