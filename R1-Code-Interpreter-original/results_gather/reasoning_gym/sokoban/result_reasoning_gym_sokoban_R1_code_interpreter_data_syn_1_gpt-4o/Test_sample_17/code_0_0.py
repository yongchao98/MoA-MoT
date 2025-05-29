from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (6, 2),  # Player position
    ((3, 2), (5, 1), (6, 3))  # Box positions
)

# Define the goal positions
goals = {(3, 1), (6, 2)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
    (1, 0), (1, 1), (1, 5),
    (2, 0), (2, 5),
    (3, 0), (3, 4), (3, 5),
    (4, 0), (4, 4), (4, 5),
    (5, 0), (5, 4), (5, 5),
    (6, 0), (6, 4), (6, 5),
    (7, 0), (7, 1), (7, 2), (7, 3), (7, 4), (7, 5)
}

# Define the possible moves
moves = {
    'L': (0, -1),
    'R': (0, 1),
    'U': (-1, 0),
    'D': (1, 0)
}

# Function to check if a position is free (not a wall or a box)
def is_free(pos, boxes):
    return pos not in walls and pos not in boxes

# Function to apply a move
def apply_move(state, move):
    player, boxes = state
    dx, dy = moves[move]
    new_player = (player[0] + dx, player[1] + dy)
    
    if new_player in boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if is_free(new_box, boxes):
            new_boxes = tuple(b if b != new_player else new_box for b in boxes)
            return new_player, new_boxes
    elif is_free(new_player, boxes):
        return new_player, boxes
    
    return None

# Function to check if the puzzle is solved
def is_solved(boxes):
    return all(box in goals for box in boxes)

# Breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        state, path = queue.popleft()
        
        if is_solved(state[1]):
            return path
        
        for move in moves:
            new_state = apply_move(state, move)
            if new_state and new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + move))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)