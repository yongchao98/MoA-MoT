from heapq import heappop, heappush

# Define the initial state of the puzzle
initial_state = {
    'player': (1, 5),
    'boxes': {(2, 2), (3, 1), (3, 4), (5, 2), (5, 5)},
    'goals': {(1, 6), (3, 2), (4, 1), (4, 4)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
              (1, 0), (1, 7), (2, 0), (2, 7), (3, 0), (3, 7), (4, 0), (4, 7),
              (5, 0), (5, 7), (6, 0), (6, 7), (7, 0), (7, 1), (7, 2), (7, 3),
              (7, 4), (7, 5), (7, 6), (7, 7)}
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Function to check if a position is within bounds and not a wall
def is_free(pos, state):
    return pos not in state['walls'] and pos not in state['boxes']

# Function to apply a move
def apply_move(state, move):
    player = state['player']
    boxes = state['boxes']
    dx, dy = moves[move]
    new_player = (player[0] + dx, player[1] + dy)
    
    if new_player in boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if is_free(new_box, state):
            new_boxes = boxes - {new_player} | {new_box}
            return {'player': new_player, 'boxes': new_boxes, 'goals': state['goals'], 'walls': state['walls']}
    elif is_free(new_player, state):
        return {'player': new_player, 'boxes': boxes, 'goals': state['goals'], 'walls': state['walls']}
    return None

# Function to check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Heuristic function: sum of Manhattan distances from boxes to nearest goals
def heuristic(state):
    total_distance = 0
    for box in state['boxes']:
        min_distance = min(abs(box[0] - goal[0]) + abs(box[1] - goal[1]) for goal in state['goals'])
        total_distance += min_distance
    return total_distance

# A* search to find the solution
def solve_sokoban(initial_state):
    queue = []
    heappush(queue, (0 + heuristic(initial_state), 0, initial_state, ""))
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while queue:
        priority, cost, current_state, path = heappop(queue)
        
        if is_solved(current_state):
            return path
        
        for move in moves:
            new_state = apply_move(current_state, move)
            if new_state:
                state_id = (new_state['player'], frozenset(new_state['boxes']))
                if state_id not in visited:
                    visited.add(state_id)
                    new_cost = cost + 1
                    new_priority = new_cost + heuristic(new_state)
                    heappush(queue, (new_priority, new_cost, new_state, path + move))
    return None

# Solve the puzzle
solution = solve_sokoban(initial_state)
print(solution)