from heapq import heappop, heappush

# Define the initial state of the puzzle
initial_state = {
    'player': (8, 9),
    'boxes': {(2, 3), (2, 4), (3, 3), (6, 3), (8, 5)},
    'goals': {(1, 2), (1, 5), (3, 2), (3, 5), (8, 2)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
              (1, 0), (1, 6), (2, 0), (2, 6), (3, 0), (3, 6),
              (4, 0), (4, 5), (4, 6), (5, 0), (5, 5), (5, 6),
              (6, 0), (6, 5), (6, 6), (7, 0), (7, 5), (7, 6),
              (8, 0), (8, 6), (9, 0), (9, 1), (9, 2), (9, 3), (9, 4), (9, 5), (9, 6)},
    'width': 7,
    'height': 10
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

def is_goal_state(state):
    return state['boxes'] == state['goals']

def move_player(state, direction):
    dx, dy = moves[direction]
    px, py = state['player']
    new_player_pos = (px + dx, py + dy)
    
    if new_player_pos in state['walls']:
        return None
    
    new_boxes = set(state['boxes'])
    if new_player_pos in new_boxes:
        new_box_pos = (new_player_pos[0] + dx, new_player_pos[1] + dy)
        if new_box_pos in state['walls'] or new_box_pos in new_boxes:
            return None
        new_boxes.remove(new_player_pos)
        new_boxes.add(new_box_pos)
    
    return {
        'player': new_player_pos,
        'boxes': new_boxes,
        'goals': state['goals'],
        'walls': state['walls'],
        'width': state['width'],
        'height': state['height']
    }

def heuristic(state):
    # Calculate the sum of Manhattan distances from each box to the nearest goal
    total_distance = 0
    for box in state['boxes']:
        min_distance = min(abs(box[0] - goal[0]) + abs(box[1] - goal[1]) for goal in state['goals'])
        total_distance += min_distance
    return total_distance

def solve_sokoban(initial_state):
    open_set = []
    heappush(open_set, (0, "", initial_state['player'], frozenset(initial_state['boxes'])))
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while open_set:
        _, path, player_pos, boxes = heappop(open_set)
        
        current_state = {
            'player': player_pos,
            'boxes': set(boxes),
            'goals': initial_state['goals'],
            'walls': initial_state['walls'],
            'width': initial_state['width'],
            'height': initial_state['height']
        }
        
        if is_goal_state(current_state):
            return path
        
        for direction in moves:
            new_state = move_player(current_state, direction)
            if new_state is None:
                continue
            
            state_signature = (new_state['player'], frozenset(new_state['boxes']))
            if state_signature not in visited:
                visited.add(state_signature)
                cost = len(path) + 1 + heuristic(new_state)
                heappush(open_set, (cost, path + direction, new_state['player'], frozenset(new_state['boxes'])))
    
    return None

solution = solve_sokoban(initial_state)
print(solution)