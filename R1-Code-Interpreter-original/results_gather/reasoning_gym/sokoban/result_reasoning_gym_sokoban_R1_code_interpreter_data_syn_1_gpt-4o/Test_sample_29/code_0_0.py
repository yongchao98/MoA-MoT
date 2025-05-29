from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 2), (2, 2), (2, 4), (2, 5),  # Box positions
    (5, 2)  # Player position
)

# Define the goal positions
goal_positions = {(1, 4), (1, 5), (2, 3), (3, 3)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8),
    (1, 0), (1, 7), (1, 8),
    (2, 0), (2, 7), (2, 8),
    (3, 0), (3, 7), (3, 8),
    (4, 0), (4, 7), (4, 8),
    (5, 0), (5, 7), (5, 8),
    (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 7), (6, 8)
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

def is_goal_state(state):
    boxes = state[:4]
    return all(box in goal_positions for box in boxes)

def is_valid_position(position, boxes):
    return position not in walls and position not in boxes

def move(state, direction):
    boxes = list(state[:4])
    player = state[4]
    dx, dy = moves[direction]
    new_player = (player[0] + dx, player[1] + dy)

    if new_player in boxes:
        # If the player is moving into a box, try to push the box
        box_index = boxes.index(new_player)
        new_box_position = (new_player[0] + dx, new_player[1] + dy)
        if is_valid_position(new_box_position, boxes):
            boxes[box_index] = new_box_position
            return tuple(boxes) + (new_player,)
    elif is_valid_position(new_player, boxes):
        # If the player is moving into an empty space
        return tuple(boxes) + (new_player,)

    return None

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