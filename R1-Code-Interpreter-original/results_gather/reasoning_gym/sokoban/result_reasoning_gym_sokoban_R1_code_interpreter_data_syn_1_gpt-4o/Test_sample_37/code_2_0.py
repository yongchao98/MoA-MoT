from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 1),  # Player position
    ((2, 2), (3, 2), (4, 2)),  # Box positions
)

# Define the goal positions
goals = {(3, 3), (3, 4), (4, 4)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
    (1, 0), (1, 5),
    (2, 0), (2, 5),
    (3, 0), (3, 5),
    (4, 0), (4, 5),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5),
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

def is_goal_state(boxes):
    return all(box in goals for box in boxes)

def is_valid_position(position, boxes):
    return position not in walls and position not in boxes

def move_player(player, direction):
    return (player[0] + direction[0], player[1] + direction[1])

def move_box(box, direction):
    return (box[0] + direction[0], box[1] + direction[1])

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()

    while queue:
        (player, boxes), path = queue.popleft()

        if is_goal_state(boxes):
            return path

        for move, direction in moves.items():
            new_player = move_player(player, direction)

            if is_valid_position(new_player, boxes):
                new_boxes = list(boxes)
                if new_player in boxes:
                    new_box = move_box(new_player, direction)
                    if is_valid_position(new_box, boxes):
                        new_boxes.remove(new_player)
                        new_boxes.append(new_box)
                    else:
                        continue

                new_state = (new_player, tuple(sorted(new_boxes)))
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move))

    return "No solution"

solution = solve_sokoban(initial_state)
print(solution)