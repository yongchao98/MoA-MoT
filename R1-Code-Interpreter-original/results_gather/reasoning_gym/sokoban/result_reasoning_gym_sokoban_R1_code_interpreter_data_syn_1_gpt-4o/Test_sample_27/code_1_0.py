from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 2),  # Box 1 position
    (2, 4),  # Box 2 position
    (3, 4),  # Box 3 position
    (3, 5)   # Player position
)

# Define the goal positions
goal_positions = {(1, 5), (2, 3), (3, 2)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
    (1, 0), (1, 6),
    (2, 0), (2, 6),
    (3, 0), (3, 6),
    (4, 0), (4, 6),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6)
}

# Define the possible moves
moves = {
    'L': (0, -1),
    'R': (0, 1),
    'U': (-1, 0),
    'D': (1, 0)
}

def is_goal_state(state):
    boxes = state[:3]
    return all(box in goal_positions for box in boxes)

def is_valid_position(position, boxes):
    return position not in walls and position not in boxes

def get_new_position(position, move):
    return (position[0] + move[0], position[1] + move[1])

def solve_sokoban():
    queue = deque([(initial_state, "")])
    visited = set()

    while queue:
        current_state, path = queue.popleft()
        boxes, player = current_state[:3], current_state[3]

        if is_goal_state(current_state):
            return path

        if current_state in visited:
            continue

        visited.add(current_state)

        for move_name, move in moves.items():
            new_player = get_new_position(player, move)

            if not is_valid_position(new_player, boxes):
                continue

            if new_player in boxes:
                box_index = boxes.index(new_player)
                new_box_position = get_new_position(new_player, move)

                if not is_valid_position(new_box_position, boxes):
                    continue

                new_boxes = list(boxes)
                new_boxes[box_index] = new_box_position
                new_state = tuple(new_boxes) + (new_player,)
            else:
                new_state = boxes + (new_player,)

            if new_state not in visited:
                queue.append((new_state, path + move_name))

    return None

solution = solve_sokoban()
print(solution)