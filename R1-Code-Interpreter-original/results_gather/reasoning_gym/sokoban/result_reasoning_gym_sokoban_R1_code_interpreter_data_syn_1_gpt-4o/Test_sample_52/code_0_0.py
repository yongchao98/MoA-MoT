from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (2, 3),  # Box 1 position
    (3, 4),  # Box 2 position
    (4, 4),  # Box 3 position
    (4, 5)   # Player position
)

# Define the goal positions
goals = {(3, 2), (4, 4)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
    (1, 0), (1, 7),
    (2, 0), (2, 7),
    (3, 0), (3, 7),
    (4, 0), (4, 7),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7)
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
    return all(box in goals for box in boxes)

def is_valid_position(position, boxes):
    return position not in walls and position not in boxes

def get_new_state(state, move):
    boxes = list(state[:3])
    player = state[3]
    dx, dy = moves[move]
    new_player = (player[0] + dx, player[1] + dy)

    if new_player in boxes:
        box_index = boxes.index(new_player)
        new_box = (boxes[box_index][0] + dx, boxes[box_index][1] + dy)
        if is_valid_position(new_box, boxes):
            boxes[box_index] = new_box
        else:
            return None
    elif not is_valid_position(new_player, boxes):
        return None

    return tuple(boxes) + (new_player,)

def solve_sokoban():
    queue = deque([(initial_state, "")])
    visited = set()

    while queue:
        current_state, path = queue.popleft()

        if is_goal_state(current_state):
            return path

        if current_state in visited:
            continue

        visited.add(current_state)

        for move in moves:
            new_state = get_new_state(current_state, move)
            if new_state and new_state not in visited:
                queue.append((new_state, path + move))

    return "No solution"

solution = solve_sokoban()
print(solution)