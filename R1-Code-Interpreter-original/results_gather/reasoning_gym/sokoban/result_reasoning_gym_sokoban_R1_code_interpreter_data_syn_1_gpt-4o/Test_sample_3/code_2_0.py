from heapq import heappop, heappush

# Define the initial state of the puzzle
initial_state = (
    (6, 5),  # Player position
    frozenset({(1, 8), (3, 7), (4, 6), (5, 5), (6, 6)})  # Box positions
)

# Define the goal positions
goals = {(2, 3), (1, 4), (2, 6), (3, 7)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9),
    (1, 0), (1, 1), (1, 7), (1, 9),
    (2, 0), (2, 7), (2, 8), (2, 9),
    (3, 0), (3, 8), (3, 9),
    (4, 0), (4, 7), (4, 8), (4, 9),
    (5, 0), (5, 6), (5, 7), (5, 8), (5, 9),
    (6, 0), (6, 7), (6, 8), (6, 9),
    (7, 0), (7, 1), (7, 2), (7, 7), (7, 8), (7, 9)
}

# Define the possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Function to check if a position is free (not a wall or a box)
def is_free(position, boxes):
    return position not in walls and position not in boxes

# Function to apply a move
def apply_move(state, move):
    player, boxes = state
    dx, dy = moves[move]
    new_player = (player[0] + dx, player[1] + dy)
    
    if new_player in boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if is_free(new_box, boxes):
            new_boxes = frozenset(b if b != new_player else new_box for b in boxes)
            return new_player, new_boxes
    elif is_free(new_player, boxes):
        return new_player, boxes
    return None

# Function to check if the puzzle is solved
def is_solved(boxes):
    return all(box in goals for box in boxes)

# Heuristic function: sum of Manhattan distances from each box to the nearest goal
def heuristic(boxes):
    return sum(min(abs(box[0] - goal[0]) + abs(box[1] - goal[1]) for goal in goals) for box in boxes)

# Function to detect deadlocks
def is_deadlock(boxes):
    for box in boxes:
        if box not in goals:
            # Check for corner deadlocks
            if ((box[0] - 1, box[1]) in walls or (box[0] + 1, box[1]) in walls) and \
               ((box[0], box[1] - 1) in walls or (box[0], box[1] + 1) in walls):
                return True
    return False

# A* search to find the solution
def solve_sokoban(initial_state):
    queue = []
    heappush(queue, (0 + heuristic(initial_state[1]), 0, initial_state, ""))
    visited = set()
    visited.add(initial_state)
    
    while queue:
        _, cost, state, path = heappop(queue)
        
        if is_solved(state[1]):
            return path
        
        for move in moves:
            new_state = apply_move(state, move)
            if new_state and new_state not in visited and not is_deadlock(new_state[1]):
                visited.add(new_state)
                new_cost = cost + 1
                heappush(queue, (new_cost + heuristic(new_state[1]), new_cost, new_state, path + move))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)