from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (5, 3),  # Player position
    ((1, 4), (2, 4), (3, 2), (4, 4)),  # Box positions
)

# Define the goal positions
goal_positions = {(1, 5), (3, 4), (4, 2)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
    (1, 0), (1, 6),
    (2, 0), (2, 6),
    (3, 0), (3, 6),
    (4, 0), (4, 6),
    (5, 0), (5, 5), (5, 6),
    (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6),
}

# Define the possible moves
moves = {
    'L': (0, -1),
    'R': (0, 1),
    'U': (-1, 0),
    'D': (1, 0),
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
    
    if new_player in walls:
        return None
    
    new_boxes = list(boxes)
    if new_player in boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if new_box in walls or new_box in boxes:
            return None
        new_boxes.remove(new_player)
        new_boxes.append(new_box)
    
    return (new_player, tuple(sorted(new_boxes)))

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal_state(current_state):
            return path
        
        for direction in moves:
            new_state = move(current_state, direction)
            if new_state and new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + direction))
    
    return None

solution = solve_sokoban(initial_state)
print(solution)