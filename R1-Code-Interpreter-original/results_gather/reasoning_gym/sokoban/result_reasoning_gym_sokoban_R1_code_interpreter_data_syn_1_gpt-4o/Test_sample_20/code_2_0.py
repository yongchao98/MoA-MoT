from collections import deque

# Define the initial state of the puzzle
initial_state = (
    ((1, 3), (1, 5), (2, 5), (3, 5), (3, 4)),  # Box positions
    (4, 3)   # Player position
)

# Define the goal positions
goal_positions = {(1, 1), (1, 2), (2, 3), (2, 4)}

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
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

def is_goal_state(box_positions):
    return all(pos in goal_positions for pos in box_positions)

def move_player(player_pos, direction):
    return (player_pos[0] + direction[0], player_pos[1] + direction[1])

def move_box(box_pos, direction):
    return (box_pos[0] + direction[0], box_pos[1] + direction[1])

def is_valid_position(pos, box_positions):
    return pos not in walls and pos not in box_positions

def solve_sokoban():
    queue = deque([(initial_state, "")])
    visited = set()

    while queue:
        state, path = queue.popleft()
        box_positions, player_pos = state

        if is_goal_state(box_positions):
            return path

        if state in visited:
            continue

        visited.add(state)

        for move, direction in moves.items():
            new_player_pos = move_player(player_pos, direction)

            if not is_valid_position(new_player_pos, box_positions):
                continue

            new_box_positions = list(box_positions)
            if new_player_pos in box_positions:
                new_box_pos = move_box(new_player_pos, direction)
                if not is_valid_position(new_box_pos, box_positions):
                    continue
                new_box_positions.remove(new_player_pos)
                new_box_positions.append(new_box_pos)

            new_state = (tuple(sorted(new_box_positions)), new_player_pos)
            if new_state not in visited:
                queue.append((new_state, path + move))

    return None

solution = solve_sokoban()
print(solution)