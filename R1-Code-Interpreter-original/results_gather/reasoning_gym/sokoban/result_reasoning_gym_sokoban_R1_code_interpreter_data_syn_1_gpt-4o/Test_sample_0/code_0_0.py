from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 1),  # Player position
    ((1, 2), (2, 1), (2, 2), (3, 2)),  # Box positions
)

# Define the goal positions
goal_positions = {(1, 3), (2, 4), (5, 3)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
    (1, 0), (1, 4), (1, 5),
    (2, 0), (2, 5),
    (3, 0), (3, 4), (3, 5),
    (4, 0), (4, 5),
    (5, 0), (5, 5),
    (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5),
}

# Define possible moves
moves = {
    'L': (0, -1),
    'R': (0, 1),
    'U': (-1, 0),
    'D': (1, 0),
}

def is_goal_state(state):
    _, boxes = state
    return all(box in goal_positions for box in boxes)

def is_valid_position(position, boxes):
    return position not in walls and position not in boxes

def move_player(state, direction):
    player, boxes = state
    move = moves[direction]
    new_player = (player[0] + move[0], player[1] + move[1])

    if not is_valid_position(new_player, boxes):
        return None

    if new_player in boxes:
        new_box = (new_player[0] + move[0], new_player[1] + move[1])
        if not is_valid_position(new_box, boxes):
            return None
        new_boxes = tuple(new_box if box == new_player else box for box in boxes)
    else:
        new_boxes = boxes

    return (new_player, new_boxes)

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()

    while queue:
        state, path = queue.popleft()

        if is_goal_state(state):
            return path

        if state in visited:
            continue

        visited.add(state)

        for direction in moves:
            new_state = move_player(state, direction)
            if new_state and new_state not in visited:
                queue.append((new_state, path + direction))

    return None

solution = solve_sokoban(initial_state)
print(solution)