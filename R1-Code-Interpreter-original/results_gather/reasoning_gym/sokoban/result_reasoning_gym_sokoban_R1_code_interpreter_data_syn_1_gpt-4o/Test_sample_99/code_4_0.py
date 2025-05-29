from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 6),  # Player position
    frozenset([(1, 3), (1, 5), (2, 2), (2, 5), (2, 6), (4, 3)]),  # Box positions
)

# Define the goal positions
goals = frozenset([(1, 4), (1, 6), (5, 5)])

# Define the walls
walls = frozenset([
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
    (1, 0), (1, 1), (1, 2), (1, 7),
    (2, 0), (2, 1), (2, 3), (2, 4), (2, 7),
    (3, 0), (3, 1), (3, 2), (3, 3), (3, 7),
    (4, 0), (4, 1), (4, 2), (4, 7),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 6), (5, 7),
    (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 7),
    (7, 0), (7, 1), (7, 2), (7, 3), (7, 4), (7, 5), (7, 6), (7, 7),
])

# Define the possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

def is_solved(boxes):
    return boxes == goals

def is_valid_position(position, boxes):
    return position not in walls and position not in boxes

def move_player(player, direction):
    return (player[0] + direction[0], player[1] + direction[1])

def move_box(box, direction):
    return (box[0] + direction[0], box[1] + direction[1])

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set([initial_state])

    while queue:
        (player, boxes), path = queue.popleft()

        if is_solved(boxes):
            return path

        for move, direction in moves.items():
            new_player = move_player(player, direction)

            if is_valid_position(new_player, boxes):
                new_state = (new_player, boxes)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move))
            elif new_player in boxes:
                new_box = move_box(new_player, direction)
                if is_valid_position(new_box, boxes):
                    new_boxes = frozenset((new_box if b == new_player else b) for b in boxes)
                    new_state = (new_player, new_boxes)
                    if new_state not in visited:
                        visited.add(new_state)
                        queue.append((new_state, path + move))

    return None

solution = solve_sokoban(initial_state)
print(solution)