import json
from heapq import heappush, heappop
from copy import deepcopy

initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}
goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

def get_neighbors(pos):
    row, col = int(pos[1]), int(pos[3])
    neighbors = []
    if row > 1: neighbors.append(f"C{row-1},{col}")
    if row < 4: neighbors.append(f"C{row+1},{col}")
    if col > 1: neighbors.append(f"C{row},{col-1}")
    if col < 8: neighbors.append(f"C{row},{col+1}")
    return neighbors

def state_hash(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def get_moves(state):
    moves = []
    # Try to move boxes that are not at their goals
    for box in state:
        if state[box] == goal_state[box]:
            continue
        curr_pos = state[box]
        for next_pos in get_neighbors(curr_pos):
            if next_pos not in state.values():
                new_state = state.copy()
                new_state[box] = next_pos
                moves.append(new_state)
    return moves

def manhattan_dist(pos1, pos2):
    r1, c1 = int(pos1[1]), int(pos1[3])
    r2, c2 = int(pos2[1]), int(pos2[3])
    return abs(r1 - r2) + abs(c1 - c2)

def heuristic(state):
    return sum(manhattan_dist(state[box], goal_state[box]) for box in state)

def find_solution():
    visited = set()
    start_h = heuristic(initial_state)
    queue = [(start_h, 0, initial_state, [initial_state])]
    visited.add(state_hash(initial_state))
    
    while queue:
        _, g, curr_state, path = heappop(queue)
        
        if curr_state == goal_state:
            return path
            
        if g > 25:  # Limit path length
            continue
            
        for next_state in get_moves(curr_state):
            state_str = state_hash(next_state)
            if state_str not in visited:
                visited.add(state_str)
                h = heuristic(next_state)
                new_path = path + [next_state]
                heappush(queue, (g + 1 + h, g + 1, next_state, new_path))
    
    return None

# Find and format solution
solution = find_solution()
if solution:
    # Convert to list of sorted dictionaries
    formatted = []
    for state in solution:
        formatted.append(dict(sorted(state.items())))
    # Print the first state to verify
    print(json.dumps(formatted))
else:
    print("No solution found")