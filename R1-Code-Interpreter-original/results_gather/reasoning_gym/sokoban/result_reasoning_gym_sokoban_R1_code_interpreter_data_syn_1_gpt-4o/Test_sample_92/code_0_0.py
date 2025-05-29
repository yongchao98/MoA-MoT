from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 4),  # Box 1 position
    (2, 4),  # Box 2 position
    (3, 2),  # Box 3 position
    (2, 1)   # Player position
)

# Define the goal positions
goal_positions = {(1, 3), (3, 3), (3, 6)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
    (1, 0), (1, 1), (1, 2), (1, 7),
    (2, 0), (2, 7),
    (3, 0), (3, 7),
    (4, 0), (4, 7),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7)
}

# Define the possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Check if a position is within bounds and not a wall
def is_valid_position(pos):
    return pos not in walls

# Check if the puzzle is solved
def is_solved(state):
    boxes = state[:3]
    return all(box in goal_positions for box in boxes)

# Perform a move
def move(state, direction):
    boxes = list(state[:3])
    player = state[3]
    move_offset = moves[direction]
    new_player = (player[0] + move_offset[0], player[1] + move_offset[1])

    if not is_valid_position(new_player):
        return None

    if new_player in boxes:
        box_index = boxes.index(new_player)
        new_box = (boxes[box_index][0] + move_offset[0], boxes[box_index][1] + move_offset[1])
        if not is_valid_position(new_box) or new_box in boxes:
            return None
        boxes[box_index] = new_box

    return tuple(boxes) + (new_player,)

# Breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, path = queue.popleft()

        if is_solved(current_state):
            return path

        for direction in moves:
            new_state = move(current_state, direction)
            if new_state and new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + direction))

    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)