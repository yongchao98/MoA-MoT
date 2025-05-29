from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (3, 2),  # Player position
    frozenset(((3, 4), (4, 3), (6, 2), (6, 4))),  # Box positions
)

# Define the goal positions
goal_positions = frozenset({(5, 2), (6, 3), (4, 5), (6, 4)})

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
    (1, 0), (1, 6),
    (2, 0), (2, 6),
    (3, 0), (3, 6),
    (4, 0), (4, 5), (4, 6),
    (5, 0), (5, 6),
    (6, 0), (6, 1), (6, 5), (6, 6),
    (7, 0), (7, 1), (7, 2), (7, 3), (7, 4), (7, 5), (7, 6),
}

# Define the possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

def is_valid_position(position, boxes):
    return position not in walls and position not in boxes

def move_player(state, move):
    player, boxes = state
    dx, dy = moves[move]
    new_player = (player[0] + dx, player[1] + dy)
    
    if not is_valid_position(new_player, boxes):
        return None
    
    if new_player in boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if not is_valid_position(new_box, boxes):
            return None
        new_boxes = frozenset(b if b != new_player else new_box for b in boxes)
    else:
        new_boxes = boxes
    
    return (new_player, new_boxes)

def is_goal_state(state):
    _, boxes = state
    return boxes == goal_positions

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    
    while queue:
        state, path = queue.popleft()
        
        if is_goal_state(state):
            return path
        
        if state in visited:
            continue
        
        visited.add(state)
        
        for move in moves:
            new_state = move_player(state, move)
            if new_state and new_state not in visited:
                queue.append((new_state, path + move))
    
    return None

solution = solve_sokoban(initial_state)
print(solution)