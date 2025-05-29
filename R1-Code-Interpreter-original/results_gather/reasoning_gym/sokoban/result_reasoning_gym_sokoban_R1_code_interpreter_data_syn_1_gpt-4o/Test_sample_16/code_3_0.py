from collections import deque

# Define the initial state of the puzzle
initial_state = (
    ((1, 4), (2, 2), (5, 3), (6, 3), (6, 5)),  # Box positions
    (6, 4)   # Player position
)

# Define the goal positions
goal_positions = {(1, 2), (2, 3), (4, 3), (6, 5)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9),
    (1, 0), (1, 7), (1, 8), (1, 9),
    (2, 0), (2, 7), (2, 8), (2, 9),
    (3, 0), (3, 7), (3, 8), (3, 9),
    (4, 0), (4, 5), (4, 7), (4, 8), (4, 9),
    (5, 0), (5, 5), (5, 7), (5, 8), (5, 9),
    (6, 0), (6, 7), (6, 8), (6, 9),
    (7, 0), (7, 1), (7, 2), (7, 3), (7, 4), (7, 5), (7, 6), (7, 7), (7, 8), (7, 9)
}

# Define the possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

def is_goal(state):
    # Check if all boxes are on goal positions
    boxes, _ = state
    return all(box in goal_positions for box in boxes)

def is_valid(state, move):
    # Check if the move is valid
    boxes, player = state
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
    # Apply the move to the state
    boxes, player = list(state[0]), state[1]
    dx, dy = moves[move]
    new_player = (player[0] + dx, player[1] + dy)
    
    if new_player in boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        boxes[boxes.index(new_player)] = new_box
    
    return (tuple(boxes), new_player)

def solve_sokoban(initial_state):
    # Solve the Sokoban puzzle using BFS
    queue = deque([(initial_state, "")])
    visited = set()
    
    while queue:
        state, path = queue.popleft()
        
        if is_goal(state):
            return path
        
        for move in moves:
            if is_valid(state, move):
                new_state = apply_move(state, move)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move))
    
    return None

solution = solve_sokoban(initial_state)
print(solution)