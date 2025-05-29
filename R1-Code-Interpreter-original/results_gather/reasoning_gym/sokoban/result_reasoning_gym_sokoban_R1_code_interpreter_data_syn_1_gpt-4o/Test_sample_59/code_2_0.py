from heapq import heappop, heappush
from itertools import count

# Define the initial state
initial_state = {
    'player': (7, 2),
    'boxes': {(1, 3), (2, 5), (4, 3), (4, 4), (5, 3), (6, 2)},
    'goals': {(1, 2), (1, 5), (1, 6), (3, 5), (5, 6), (6, 1)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
              (1, 0), (1, 7), (2, 0), (2, 7), (3, 0), (3, 7), (4, 0), (4, 6), (4, 7),
              (5, 0), (5, 7), (6, 0), (6, 7), (7, 0), (7, 7), (8, 0), (8, 1), (8, 2),
              (8, 3), (8, 4), (8, 5), (8, 6), (8, 7)},
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
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
        'walls': state['walls'],
    }

# Check if the state is a goal state
def is_goal_state(state):
    return state['boxes'] == state['goals']

# Heuristic function: sum of Manhattan distances from boxes to nearest goals
def heuristic(state):
    total_distance = 0
    for box in state['boxes']:
        min_distance = min(abs(box[0] - goal[0]) + abs(box[1] - goal[1]) for goal in state['goals'])
        total_distance += min_distance
    return total_distance

# Perform A* search to find the solution
def solve_sokoban(initial_state):
    open_set = []
    counter = count()
    heappush(open_set, (0 + heuristic(initial_state), next(counter), 0, initial_state, ""))
    visited = set()
    
    while open_set:
        _, _, cost, current_state, path = heappop(open_set)
        
        if is_goal_state(current_state):
            return path
        
        state_key = (current_state['player'], frozenset(current_state['boxes']))
        if state_key in visited:
            continue
        visited.add(state_key)
        
        for move_name, move in moves.items():
            if is_valid_move(current_state, move):
                new_state = apply_move(current_state, move)
                new_cost = cost + 1
                heappush(open_set, (new_cost + heuristic(new_state), next(counter), new_cost, new_state, path + move_name))
    
    return None

# Solve the puzzle
solution = solve_sokoban(initial_state)
print(solution)