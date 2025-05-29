from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (1, 2),  # Box on goal
    (3, 4),  # Goal
    (5, 2),  # Goal
    (4, 4),  # Box
    (4, 5),  # Box
    (5, 5),  # Box
    (4, 6)   # Player
)

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
    (1, 0), (1, 1), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7),
    (2, 0), (2, 1), (2, 2), (2, 7),
    (3, 0), (3, 7),
    (4, 0), (4, 7),
    (5, 0), (5, 7),
    (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 7)
}

# Define the directions
directions = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Check if a position is within bounds and not a wall
def is_valid(pos):
    return pos not in walls

# Check if the puzzle is solved
def is_solved(state):
    goals = state[:2]
    boxes = state[3:6]
    return all(box in goals for box in boxes)

# Perform a move
def move(state, direction):
    player = state[6]
    dx, dy = directions[direction]
    new_player = (player[0] + dx, player[1] + dy)
    
    if not is_valid(new_player):
        return None
    
    new_boxes = list(state[3:6])
    if new_player in new_boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if not is_valid(new_box) or new_box in new_boxes:
            return None
        new_boxes[new_boxes.index(new_player)] = new_box
    
    return state[:3] + tuple(new_boxes) + (new_player,)

# Solve the puzzle using BFS
def solve(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    
    while queue:
        state, path = queue.popleft()
        
        if is_solved(state):
            return path
        
        for direction in directions:
            new_state = move(state, direction)
            if new_state and new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + direction))
    
    return None

# Find the solution
solution = solve(initial_state)
print(solution)