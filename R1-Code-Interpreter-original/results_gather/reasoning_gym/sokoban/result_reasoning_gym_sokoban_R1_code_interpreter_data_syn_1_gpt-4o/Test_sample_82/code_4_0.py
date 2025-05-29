from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (5, 3),  # Box
    (6, 2),  # Box
    (8, 3),  # Box
    (8, 2),  # Player
)

# Define the goal positions
goals = {(2, 1), (6, 3), (8, 4)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
    (1, 0), (1, 1), (1, 3), (1, 4), (1, 5), (1, 6),
    (2, 0), (2, 2), (2, 3), (2, 4), (2, 6),
    (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6),
    (4, 0), (4, 1), (4, 4), (4, 5), (4, 6),
    (5, 0), (5, 1), (5, 4), (5, 6),
    (6, 0), (6, 1), (6, 4), (6, 5), (6, 6),
    (7, 0), (7, 4), (7, 5), (7, 6),
    (8, 0), (8, 5), (8, 6),
    (9, 0), (9, 1), (9, 2), (9, 3), (9, 4), (9, 5), (9, 6),
}

# Define the possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

# Function to check if a position is within bounds and not a wall
def is_valid(pos):
    return pos not in walls

# Function to apply a move
def apply_move(state, move):
    player = state[-1]
    new_player = (player[0] + move[0], player[1] + move[1])
    if not is_valid(new_player):
        return None
    new_boxes = list(state[:-1])
    if new_player in new_boxes:
        new_box = (new_player[0] + move[0], new_player[1] + move[1])
        if not is_valid(new_box) or new_box in new_boxes:
            return None
        new_boxes[new_boxes.index(new_player)] = new_box
    return tuple(new_boxes + [new_player])

# Function to check if the puzzle is solved
def is_solved(state):
    return all(box in goals for box in state[:-1])

# Breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    while queue:
        state, path = queue.popleft()
        if is_solved(state):
            return path
        if state in visited:
            continue
        visited.add(state)
        for move_name, move in moves.items():
            new_state = apply_move(state, move)
            if new_state and new_state not in visited:
                queue.append((new_state, path + move_name))
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)