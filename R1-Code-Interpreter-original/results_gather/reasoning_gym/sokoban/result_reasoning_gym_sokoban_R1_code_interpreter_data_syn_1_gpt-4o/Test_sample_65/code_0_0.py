from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 4),  # Box positions
    (3, 3),
    (3, 4),
    (5, 3),
    (6, 2),
    (6, 4),
    (7, 5)   # Player position
)

# Define the goal positions
goal_positions = {(1, 3), (3, 2), (4, 2), (5, 2), (7, 2)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
    (1, 0), (1, 1), (1, 2), (1, 5), (1, 6),
    (2, 0), (2, 1), (2, 6),
    (3, 0), (3, 6),
    (4, 0), (4, 6),
    (5, 0), (5, 6),
    (6, 0), (6, 6),
    (7, 0), (7, 6),
    (8, 0), (8, 1), (8, 2), (8, 3), (8, 4), (8, 5), (8, 6)
}

# Define the possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Function to check if a position is within bounds and not a wall
def is_valid_position(pos):
    return pos not in walls

# Function to apply a move to a position
def apply_move(pos, move):
    return (pos[0] + move[0], pos[1] + move[1])

# Function to check if the puzzle is solved
def is_solved(state):
    box_positions = state[:-1]
    return all(pos in goal_positions for pos in box_positions)

# Breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, path = queue.popleft()
        box_positions = current_state[:-1]
        player_position = current_state[-1]

        if is_solved(current_state):
            return path

        for move, delta in moves.items():
            new_player_position = apply_move(player_position, delta)

            if not is_valid_position(new_player_position):
                continue

            if new_player_position in box_positions:
                new_box_position = apply_move(new_player_position, delta)

                if not is_valid_position(new_box_position) or new_box_position in box_positions:
                    continue

                new_box_positions = tuple(
                    new_box_position if pos == new_player_position else pos
                    for pos in box_positions
                )
            else:
                new_box_positions = box_positions

            new_state = new_box_positions + (new_player_position,)

            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + move))

    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)