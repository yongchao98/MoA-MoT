from collections import deque

# Define the initial state
initial_state = (
    (1, 6),  # Player position
    ((1, 5), (2, 6), (3, 1))  # Box positions
)

# Define the goal positions
goal_positions = {(1, 1), (1, 4), (4, 6)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
    (1, 0), (1, 2), (1, 7),
    (2, 0), (2, 7),
    (3, 0), (3, 7),
    (4, 0), (4, 7),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7)
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

def is_goal_state(state):
    _, boxes = state
    return all(box in goal_positions for box in boxes)

def is_valid_position(position, boxes):
    return position not in walls and position not in boxes

def move(state, direction):
    player, boxes = state
    dx, dy = moves[direction]
    new_player = (player[0] + dx, player[1] + dy)
    
    if not is_valid_position(new_player, boxes):
        return None
    
    new_boxes = set(boxes)
    if new_player in boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if not is_valid_position(new_box, boxes):
            return None
        new_boxes.remove(new_player)
        new_boxes.add(new_box)
    
    return (new_player, tuple(sorted(new_boxes)))

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal_state(current_state):
            return path
        
        if current_state in visited:
            continue
        
        visited.add(current_state)
        
        for direction in moves:
            new_state = move(current_state, direction)
            if new_state and new_state not in visited:
                queue.append((new_state, path + direction))
    
    return None

solution = solve_sokoban(initial_state)
print(solution)