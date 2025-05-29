from heapq import heappop, heappush
from itertools import count

# Define the initial state of the puzzle
initial_state = {
    'player': (9, 5),
    'boxes': {(1, 1), (3, 4), (4, 4), (7, 4), (8, 4), (8, 5)},
    'goals': {(1, 5), (4, 2), (5, 4), (6, 4), (6, 5), (6, 6)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
              (1, 0), (1, 2), (1, 6),
              (2, 0), (2, 6),
              (3, 0), (3, 6),
              (4, 0), (4, 6),
              (5, 0), (5, 6),
              (6, 0), (6, 6),
              (7, 0), (7, 6),
              (8, 0), (8, 6),
              (9, 0), (9, 6),
              (10, 0), (10, 1), (10, 2), (10, 3), (10, 4), (10, 5), (10, 6)}
}

# Define possible moves
moves = {'L': (0, -1), 'R': (0, 1), 'U': (-1, 0), 'D': (1, 0)}

# Function to check if a position is free (not a wall or a box)
def is_free(position, boxes, walls):
    return position not in walls and position not in boxes

# Function to apply a move
def apply_move(state, move):
    player = state['player']
    boxes = state['boxes']
    dx, dy = moves[move]
    new_player = (player[0] + dx, player[1] + dy)
    
    if new_player in boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if is_free(new_box, boxes, state['walls']):
            new_boxes = boxes - {new_player} | {new_box}
            return {'player': new_player, 'boxes': new_boxes, 'goals': state['goals'], 'walls': state['walls']}
    elif is_free(new_player, boxes, state['walls']):
        return {'player': new_player, 'boxes': boxes, 'goals': state['goals'], 'walls': state['walls']}
    return None

# Function to check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Enhanced heuristic function
def heuristic(state):
    total_distance = 0
    for box in state['boxes']:
        min_distance = min(abs(box[0] - goal[0]) + abs(box[1] - goal[1]) for goal in state['goals'])
        total_distance += min_distance
    # Add distance from player to nearest box
    player = state['player']
    player_distance = min(abs(player[0] - box[0]) + abs(player[1] - box[1]) for box in state['boxes'])
    return total_distance + player_distance

# A* search to find the solution
def solve_sokoban(initial_state):
    open_set = []
    counter = count()  # Unique sequence count
    heappush(open_set, (0, next(counter), initial_state, ""))
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while open_set:
        _, _, current_state, path = heappop(open_set)
        
        if is_solved(current_state):
            return path
        
        for move in moves:
            new_state = apply_move(current_state, move)
            if new_state:
                state_id = (new_state['player'], frozenset(new_state['boxes']))
                if state_id not in visited:
                    visited.add(state_id)
                    cost = len(path) + 1 + heuristic(new_state)
                    heappush(open_set, (cost, next(counter), new_state, path + move))
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)