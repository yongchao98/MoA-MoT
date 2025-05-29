from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (1, 7),
    'boxes': {(3, 4), (4, 4), (6, 5)},
    'goals': {(1, 5), (3, 3), (4, 3), (4, 5)}
}

# Define the walls and the grid size
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8),
    (1, 0), (1, 2), (1, 3), (1, 8),
    (2, 0), (2, 1), (2, 2), (2, 3), (2, 8),
    (3, 0), (3, 1), (3, 2), (3, 8),
    (4, 0), (4, 1), (4, 2), (4, 8),
    (5, 0), (5, 1), (5, 2), (5, 8),
    (6, 0), (6, 1), (6, 8),
    (7, 0), (7, 8),
    (8, 0), (8, 1), (8, 2), (8, 3), (8, 4), (8, 5), (8, 6), (8, 7), (8, 8)
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Function to check if a position is within bounds and not a wall
def is_valid_position(pos):
    return pos not in walls

# Function to apply a move
def apply_move(state, move):
    player = state['player']
    boxes = state['boxes']
    dx, dy = moves[move]
    new_player = (player[0] + dx, player[1] + dy)
    
    if new_player in walls:
        return None
    
    new_boxes = set(boxes)
    if new_player in boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if new_box in walls or new_box in boxes:
            return None
        new_boxes.remove(new_player)
        new_boxes.add(new_box)
    
    return {
        'player': new_player,
        'boxes': new_boxes,
        'goals': state['goals']
    }

# Function to check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
        
        for move in moves:
            new_state = apply_move(current_state, move)
            if new_state is None:
                continue
            
            state_key = (new_state['player'], frozenset(new_state['boxes']))
            if state_key not in visited:
                visited.add(state_key)
                queue.append((new_state, path + move))
    
    return None

# Solve the puzzle
solution = solve_sokoban(initial_state)
print(solution)