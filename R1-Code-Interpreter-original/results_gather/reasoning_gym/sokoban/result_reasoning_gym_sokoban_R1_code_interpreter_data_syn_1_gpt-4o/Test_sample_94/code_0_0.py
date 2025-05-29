from collections import deque

# Define the initial state of the puzzle
initial_state = (
    (2, 4),  # Box 1 position
    (2, 6),  # Box 2 position
    (3, 5),  # Box 3 position
    (4, 4)   # Player position
)

# Define the goal positions
goals = {(2, 3), (3, 3), (4, 7)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9),
    (1, 0), (1, 1), (1, 8), (1, 9),
    (2, 0), (2, 1), (2, 8), (2, 9),
    (3, 0), (3, 1), (3, 8), (3, 9),
    (4, 0), (4, 1), (4, 8), (4, 9),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9)
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

def is_goal(state):
    boxes = state[:3]
    return all(box in goals for box in boxes)

def is_valid_move(state, move):
    boxes, player = state[:3], state[3]
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
    boxes, player = list(state[:3]), state[3]
    dx, dy = moves[move]
    new_player = (player[0] + dx, player[1] + dy)
    
    if new_player in boxes:
        box_index = boxes.index(new_player)
        boxes[box_index] = (new_player[0] + dx, new_player[1] + dy)
    
    return tuple(boxes + [new_player])

def solve_sokoban():
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal(current_state):
            return path
        
        for move in moves:
            if is_valid_move(current_state, move):
                new_state = apply_move(current_state, move)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move))
    
    return "No solution found"

solution = solve_sokoban()
print(solution)