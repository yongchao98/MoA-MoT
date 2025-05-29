import heapq

# Define the initial state of the puzzle
initial_state = {
    'player': (3, 4),
    'boxes': {(1, 2), (2, 5), (2, 7), (4, 2), (5, 5)},
    'goals': {(1, 4), (2, 2), (2, 6), (3, 8), (4, 4), (5, 4)},
    'walls': {(0, i) for i in range(10)} | {(6, i) for i in range(10)} | {(i, 0) for i in range(7)} | {(i, 9) for i in range(7)}
}

# Define the possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Function to check if a move is valid
def is_valid_move(state, move):
    player_pos = state['player']
    new_player_pos = (player_pos[0] + move[0], player_pos[1] + move[1])
    
    # Check if the new player position is a wall
    if new_player_pos in state['walls']:
        return False
    
    # Check if the new player position is a box
    if new_player_pos in state['boxes']:
        # Calculate the new box position
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        # Check if the new box position is a wall or another box
        if new_box_pos in state['walls'] or new_box_pos in state['boxes']:
            return False
    
    return True

# Function to apply a move to the state
def apply_move(state, move):
    player_pos = state['player']
    new_player_pos = (player_pos[0] + move[0], player_pos[1] + move[1])
    
    new_boxes = set(state['boxes'])
    
    # Check if the new player position is a box
    if new_player_pos in state['boxes']:
        # Calculate the new box position
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        # Move the box
        new_boxes.remove(new_player_pos)
        new_boxes.add(new_box_pos)
    
    return {
        'player': new_player_pos,
        'boxes': new_boxes,
        'goals': state['goals'],
        'walls': state['walls']
    }

# Function to check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Heuristic function: sum of Manhattan distances from each box to the nearest goal
def heuristic(state):
    total_distance = 0
    for box in state['boxes']:
        min_distance = min(abs(box[0] - goal[0]) + abs(box[1] - goal[1]) for goal in state['goals'])
        total_distance += min_distance
    return total_distance

# A* search to find the solution
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
                    heapq.heappush(priority_queue, (new_cost + heuristic(new_state), new_cost, state_key, path + move_name))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)