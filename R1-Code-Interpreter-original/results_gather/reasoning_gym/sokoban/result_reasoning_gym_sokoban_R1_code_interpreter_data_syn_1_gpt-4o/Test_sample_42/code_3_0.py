from heapq import heappop, heappush

# Define the initial state of the puzzle
initial_state = {
    'player': (1, 7),
    'boxes': {(1, 6), (2, 4), (2, 6), (2, 7), (2, 8), (3, 6)},
    'goals': {(1, 3), (2, 5), (3, 4), (3, 7), (4, 9), (5, 2)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9),
              (1, 0), (1, 9), (2, 0), (2, 9), (3, 0), (3, 9), (4, 0), (4, 4), (4, 5), (4, 6), (4, 9),
              (5, 0), (5, 3), (5, 4), (5, 5), (5, 6), (5, 9), (6, 0), (6, 4), (6, 5), (6, 6), (6, 9),
              (7, 0), (7, 5), (7, 6), (7, 7), (7, 8), (7, 9), (8, 0), (8, 5), (8, 6), (8, 7), (8, 8), (8, 9),
              (9, 0), (9, 1), (9, 2), (9, 3), (9, 4), (9, 5), (9, 6), (9, 7), (9, 8), (9, 9)}
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Check if the move is valid
def is_valid_move(state, move):
    player_pos = state['player']
    new_player_pos = (player_pos[0] + move[0], player_pos[1] + move[1])
    
    if new_player_pos in state['walls']:
        return False
    
    if new_player_pos in state['boxes']:
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        if new_box_pos in state['walls'] or new_box_pos in state['boxes']:
            return False
    
    return True

# Apply the move to the state
def apply_move(state, move):
    player_pos = state['player']
    new_player_pos = (player_pos[0] + move[0], player_pos[1] + move[1])
    
    new_boxes = set(state['boxes'])
    if new_player_pos in new_boxes:
        new_boxes.remove(new_player_pos)
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        new_boxes.add(new_box_pos)
    
    return {
        'player': new_player_pos,
        'boxes': new_boxes,
        'goals': state['goals'],
        'walls': state['walls']
    }

# Check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Calculate the heuristic cost (Manhattan distance)
def heuristic(state):
    total_distance = 0
    for box in state['boxes']:
        min_distance = min(abs(box[0] - goal[0]) + abs(box[1] - goal[1]) for goal in state['goals'])
        total_distance += min_distance
    return total_distance

# Convert state to a unique string representation
def state_to_string(state):
    return f"{state['player']}:{sorted(state['boxes'])}"

# Perform an A* search to find the solution
def solve_sokoban(initial_state):
    open_set = []
    heappush(open_set, (0, 0, initial_state, ""))
    visited = set()
    
    while open_set:
        _, _, current_state, path = heappop(open_set)
        
        if is_solved(current_state):
            return path
        
        for move_name, move in moves.items():
            if is_valid_move(current_state, move):
                new_state = apply_move(current_state, move)
                state_key = state_to_string(new_state)
                
                if state_key not in visited:
                    visited.add(state_key)
                    cost = len(path) + 1 + heuristic(new_state)
                    heappush(open_set, (cost, len(path) + 1, new_state, path + move_name))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)