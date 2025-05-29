from heapq import heappop, heappush

# Define the initial state
initial_state = {
    'player': (3, 1),
    'boxes': {(2, 1), (3, 2), (4, 1), (5, 2), (6, 7)},
    'goals': {(1, 7), (2, 5), (3, 6), (4, 6), (5, 6), (6, 6)}
}

# Define the grid
grid = [
    "+ + + + + + + + +",
    "+ - - - - - - X +",
    "+ @ - - - X $ - +",
    "+ % @ - - - @ - +",
    "+ @ - - - @ X - +",
    "+ - @ - - - X X +",
    "+ X + - - - - @ +",
    "+ - - - - - - - +",
    "+ + + + + + + + +"
]

# Convert grid to a more usable format
grid = [row.split() for row in grid]

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Function to check if a position is within bounds and not a wall
def is_valid(pos):
    x, y = pos
    return 0 <= x < len(grid) and 0 <= y < len(grid[0]) and grid[x][y] != '+'

# Function to apply a move
def apply_move(state, move):
    px, py = state['player']
    dx, dy = moves[move]
    new_player = (px + dx, py + dy)
    
    if not is_valid(new_player):
        return None
    
    new_boxes = set(state['boxes'])
    if new_player in new_boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if not is_valid(new_box) or new_box in new_boxes:
            return None
        new_boxes.remove(new_player)
        new_boxes.add(new_box)
    
    return {
        'player': new_player,
        'boxes': new_boxes,
        'goals': state['goals']
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
    open_set = []
    heappush(open_set, (0 + heuristic(initial_state), 0, initial_state, ""))
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while open_set:
        _, cost, current_state, path = heappop(open_set)
        
        if is_solved(current_state):
            return path
        
        for move in moves:
            new_state = apply_move(current_state, move)
            if new_state is None:
                continue
            
            state_id = (new_state['player'], frozenset(new_state['boxes']))
            if state_id not in visited:
                visited.add(state_id)
                new_cost = cost + 1
                heappush(open_set, (new_cost + heuristic(new_state), new_cost, new_state, path + move))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)