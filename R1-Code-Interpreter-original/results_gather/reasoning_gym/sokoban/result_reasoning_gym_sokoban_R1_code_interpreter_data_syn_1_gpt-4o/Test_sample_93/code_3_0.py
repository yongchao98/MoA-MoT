from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 2),  # Box 1
    (2, 4),  # Box 2
    (4, 1),  # Box 3
    (4, 2),  # Box 4
    (5, 2),  # Player
)

# Define the goal positions
goals = {(2, 3), (5, 3)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
    (1, 0), (1, 1), (1, 6),
    (2, 0), (2, 6),
    (3, 0), (3, 6),
    (4, 0), (4, 6),
    (5, 0), (5, 4), (5, 5), (5, 6),
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
    boxes = state[:4]
    return all(box in goals for box in boxes)

def is_valid_move(state, move):
    boxes = state[:4]
    player = state[4]
    dx, dy = moves[move]
    new_player = (player[0] + dx, player[1] + dy)
    
    if new_player in walls:
        return False
    
    if new_player in boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if new_box in walls or new_box in boxes:
            return False
    
    return True

def apply_move(state, move):
    boxes = list(state[:4])
    player = state[4]
    dx, dy = moves[move]
    new_player = (player[0] + dx, player[1] + dy)
    
    if new_player in boxes:
        box_index = boxes.index(new_player)
        boxes[box_index] = (new_player[0] + dx, new_player[1] + dy)
    
    return tuple(boxes + [new_player])

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal_state(current_state):
            return path
        
        for move in moves:
            if is_valid_move(current_state, move):
                new_state = apply_move(current_state, move)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move))
    
    return None

solution = solve_sokoban(initial_state)
print(solution)