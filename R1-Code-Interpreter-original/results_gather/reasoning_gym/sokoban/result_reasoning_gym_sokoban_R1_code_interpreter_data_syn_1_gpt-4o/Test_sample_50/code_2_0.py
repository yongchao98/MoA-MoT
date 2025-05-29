from heapq import heappop, heappush

# Define the initial state of the puzzle
initial_state = {
    'player': (7, 6),
    'boxes': {(2, 2), (3, 5), (6, 3), (6, 5)},
    'goals': {(2, 3), (5, 2), (6, 4), (7, 3)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
              (1, 0), (1, 6),
              (2, 0), (2, 6),
              (3, 0), (3, 6),
              (4, 0), (4, 6),
              (5, 0), (5, 6),
              (6, 0), (6, 6),
              (7, 0), (7, 6),
              (8, 0), (8, 1), (8, 2), (8, 3), (8, 4), (8, 5), (8, 6)}
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
    px, py = state['player']
    dx, dy = moves[direction]
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
        'walls': state['walls']
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
    heappush(open_set, (0, initial_state, ""))
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while open_set:
        _, current_state, path = heappop(open_set)
        
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
                heappush(open_set, (cost, new_state, path + direction))
    
    return None

solution = solve_sokoban(initial_state)
print(solution)