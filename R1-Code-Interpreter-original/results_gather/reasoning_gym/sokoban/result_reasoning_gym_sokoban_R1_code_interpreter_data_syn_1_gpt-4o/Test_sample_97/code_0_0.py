from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (8, 2),  # Player position
    frozenset([(2, 2), (3, 3), (3, 4), (5, 2), (7, 2), (7, 3), (7, 7)])  # Box positions
)

# Define the goal positions
goals = frozenset([(2, 3), (2, 4), (3, 5), (5, 4), (7, 5), (8, 8)])

# Define the walls
walls = frozenset([
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9),
    (1, 0), (1, 9),
    (2, 0), (2, 9),
    (3, 0), (3, 9),
    (4, 0), (4, 9),
    (5, 0), (5, 9),
    (6, 0), (6, 9),
    (7, 0), (7, 9),
    (8, 0), (8, 9),
    (9, 0), (9, 1), (9, 2), (9, 3), (9, 4), (9, 5), (9, 6), (9, 7), (9, 8), (9, 9)
])

# Define possible moves
moves = {
    'L': (0, -1),
    'R': (0, 1),
    'U': (-1, 0),
    'D': (1, 0)
}

def is_valid_move(player, box_positions, move):
    new_player = (player[0] + move[0], player[1] + move[1])
    if new_player in walls:
        return False
    if new_player in box_positions:
        new_box = (new_player[0] + move[0], new_player[1] + move[1])
        if new_box in walls or new_box in box_positions:
            return False
    return True

def apply_move(player, box_positions, move):
    new_player = (player[0] + move[0], player[1] + move[1])
    new_box_positions = set(box_positions)
    if new_player in box_positions:
        new_box = (new_player[0] + move[0], new_player[1] + move[1])
        new_box_positions.remove(new_player)
        new_box_positions.add(new_box)
    return new_player, frozenset(new_box_positions)

def is_solved(box_positions):
    return box_positions == goals

def bfs_solve(initial_state):
    queue = deque([(initial_state, "")])
    visited = set([initial_state])
    while queue:
        (player, box_positions), path = queue.popleft()
        if is_solved(box_positions):
            return path
        for move_name, move in moves.items():
            if is_valid_move(player, box_positions, move):
                new_state = apply_move(player, box_positions, move)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move_name))
    return None

solution = bfs_solve(initial_state)
print(solution)