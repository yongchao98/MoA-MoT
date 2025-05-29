from heapq import heappop, heappush

# Define the initial state
initial_state = {
    'player': (6, 3),
    'boxes': frozenset({(1, 7), (2, 6), (3, 2), (3, 4), (3, 5), (5, 2), (5, 6)}),
    'goals': {(1, 3), (2, 6), (3, 3), (3, 7), (5, 4), (6, 3), (7, 6)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9),
              (1, 0), (1, 9), (2, 0), (2, 9), (3, 0), (3, 9), (4, 0), (4, 9), (5, 0), (5, 9),
              (6, 0), (6, 9), (7, 0), (7, 9), (8, 0), (8, 9), (9, 0), (9, 1), (9, 2), (9, 3),
              (9, 4), (9, 5), (9, 6), (9, 7), (9, 8), (9, 9)},
}

# Define possible moves
moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

# Function to check if a move is valid
def is_valid_move(player, boxes, move, walls):
    new_player_pos = (player[0] + move[0], player[1] + move[1])
    
    # Check if new player position is a wall
    if new_player_pos in walls:
        return False
    
    # Check if new player position is a box
    if new_player_pos in boxes:
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        # Check if new box position is a wall or another box
        if new_box_pos in walls or new_box_pos in boxes:
            return False
    
    return True

# Function to apply a move
def apply_move(player, boxes, move):
    new_player_pos = (player[0] + move[0], player[1] + move[1])
    
    new_boxes = set(boxes)
    if new_player_pos in new_boxes:
        new_boxes.remove(new_player_pos)
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        new_boxes.add(new_box_pos)
    
    return new_player_pos, frozenset(new_boxes)

# Function to check if the puzzle is solved
def is_solved(boxes, goals):
    return boxes == goals

# Improved heuristic function
def heuristic(boxes, goals):
    total_distance = 0
    for box in boxes:
        min_distance = min(abs(box[0] - goal[0]) + abs(box[1] - goal[1]) for goal in goals)
        total_distance += min_distance
    # Add a penalty for boxes not on goals
    total_distance += len(goals) - len(boxes.intersection(goals))
    return total_distance

# A* search to find the solution
def solve_sokoban(initial_state):
    open_set = []
    player = initial_state['player']
    boxes = initial_state['boxes']
    goals = initial_state['goals']
    walls = initial_state['walls']
    
    heappush(open_set, (0, player, boxes, ""))
    visited = set()
    
    while open_set:
        _, current_player, current_boxes, path = heappop(open_set)
        
        if is_solved(current_boxes, goals):
            return path
        
        for move_key, move in moves.items():
            if is_valid_move(current_player, current_boxes, move, walls):
                new_player, new_boxes = apply_move(current_player, current_boxes, move)
                state_key = (new_player, new_boxes)
                if state_key not in visited:
                    visited.add(state_key)
                    cost = len(path) + 1 + heuristic(new_boxes, goals)
                    heappush(open_set, (cost, new_player, new_boxes, path + move_key))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)