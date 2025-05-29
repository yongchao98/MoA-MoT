from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 1),  # Player position
    frozenset([(1, 2), (2, 1), (2, 2), (2, 3), (3, 3)]),  # Box positions
)

# Define the goal positions
goal_positions = frozenset([(1, 3), (2, 3), (3, 1), (3, 2), (3, 3)])

# Define the walls
walls = frozenset([
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
    (1, 0), (1, 7),
    (2, 0), (2, 7),
    (3, 0), (3, 4), (3, 5), (3, 6), (3, 7),
    (4, 0), (4, 2), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7),
    (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 7),
    (7, 0), (7, 1), (7, 2), (7, 5), (7, 6), (7, 7),
    (8, 0), (8, 1), (8, 2), (8, 3), (8, 4), (8, 5), (8, 6), (8, 7),
])

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

def is_goal(state):
    _, boxes = state
    return boxes == goal_positions

def move_player(state, direction):
    player, boxes = state
    dx, dy = moves[direction]
    new_player = (player[0] + dx, player[1] + dy)
    
    if new_player in walls:
        return None
    
    if new_player in boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if new_box in walls or new_box in boxes:
            return None
        new_boxes = boxes - {new_player} | {new_box}
    else:
        new_boxes = boxes
    
    return (new_player, frozenset(new_boxes))

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal(current_state):
            return path
        
        for direction in moves:
            new_state = move_player(current_state, direction)
            if new_state and new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + direction))
    
    return None

solution = solve_sokoban(initial_state)
print(solution)