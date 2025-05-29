import heapq

# Define the initial state of the puzzle
initial_state = {
    'player': (3, 5),
    'boxes': {(2, 5), (3, 4), (4, 2), (8, 4)},
    'goals': {(1, 5), (4, 4), (7, 1), (8, 5)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
              (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 6),
              (2, 0), (2, 1), (2, 2), (2, 3), (2, 6),
              (3, 0), (3, 1), (3, 2), (3, 6),
              (4, 0), (4, 1), (4, 6),
              (5, 0), (5, 4), (5, 5), (5, 6),
              (6, 0), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6),
              (7, 0), (7, 2), (7, 3), (7, 4), (7, 5), (7, 6),
              (8, 0), (8, 2), (8, 3), (8, 5), (8, 6),
              (9, 0), (9, 1), (9, 2), (9, 3), (9, 4), (9, 5), (9, 6)},
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

# Heuristic function: sum of Manhattan distances from each box to the nearest goal
def heuristic(state):
    total_distance = 0
    for box in state['boxes']:
        min_distance = min(abs(box[0] - goal[0]) + abs(box[1] - goal[1]) for goal in state['goals'])
        total_distance += min_distance
    return total_distance

# Check if a move is valid
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

# Apply a move to the state
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
        'walls': state['walls'],
    }

# Check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Perform an A* search to find the solution
def solve_sokoban(initial_state):
    priority_queue = []
    heapq.heappush(priority_queue, (0 + heuristic(initial_state), 0, initial_state, ""))
    visited = set()
    
    while priority_queue:
        _, cost, current_state, path = heapq.heappop(priority_queue)
        
        if is_solved(current_state):
            return path
        
        for move_name, move in moves.items():
            if is_valid_move(current_state, move):
                new_state = apply_move(current_state, move)
                state_key = (new_state['player'], frozenset(new_state['boxes']))
                
                if state_key not in visited:
                    visited.add(state_key)
                    new_cost = cost + 1
                    heapq.heappush(priority_queue, (new_cost + heuristic(new_state), new_cost, new_state, path + move_name))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)