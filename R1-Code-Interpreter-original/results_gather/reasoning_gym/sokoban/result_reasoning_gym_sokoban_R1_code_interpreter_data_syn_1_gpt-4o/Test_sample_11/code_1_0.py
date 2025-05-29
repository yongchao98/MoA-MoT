from heapq import heappop, heappush
from itertools import count

# Define the initial state of the puzzle
initial_state = {
    'player': (1, 4),
    'boxes': {(2, 2), (2, 3), (2, 4), (2, 5)},
    'goals': {(3, 1), (5, 1), (5, 4), (5, 5)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
              (1, 0), (1, 5),
              (2, 0), (2, 5),
              (3, 0), (3, 5),
              (4, 0), (4, 5),
              (5, 0), (5, 5),
              (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5)},
    'goals_with_boxes': {(5, 3)}
}

# Define the possible moves
moves = {
    'L': (0, -1),
    'R': (0, 1),
    'U': (-1, 0),
    'D': (1, 0)
}

# Function to check if a position is free (not a wall or a box)
def is_free(position, boxes, walls):
    return position not in boxes and position not in walls

# Function to apply a move
def apply_move(state, move):
    player = state['player']
    boxes = state['boxes'].copy()
    dx, dy = moves[move]
    new_player = (player[0] + dx, player[1] + dy)
    
    if new_player in boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if is_free(new_box, boxes, state['walls']):
            boxes.remove(new_player)
            boxes.add(new_box)
            return {'player': new_player, 'boxes': boxes, 'goals': state['goals'], 'walls': state['walls']}
    elif is_free(new_player, boxes, state['walls']):
        return {'player': new_player, 'boxes': boxes, 'goals': state['goals'], 'walls': state['walls']}
    return None

# Function to check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Heuristic function: sum of Manhattan distances from boxes to nearest goals
def heuristic(state):
    return sum(min(abs(bx - gx) + abs(by - gy) for gx, gy in state['goals']) for bx, by in state['boxes'])

# A* search to find the solution
def solve_sokoban(initial_state):
    queue = []
    counter = count()
    heappush(queue, (0, next(counter), initial_state, ""))
    visited = set()
    
    while queue:
        _, _, state, path = heappop(queue)
        if is_solved(state):
            return path
        if (state['player'], frozenset(state['boxes'])) in visited:
            continue
        visited.add((state['player'], frozenset(state['boxes'])))
        
        for move in moves:
            new_state = apply_move(state, move)
            if new_state:
                new_path = path + move
                new_cost = len(new_path) + heuristic(new_state)
                heappush(queue, (new_cost, next(counter), new_state, new_path))
    return None

solution = solve_sokoban(initial_state)
print(solution)