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

# Define moves in four directions
def get_adjacent_cells(cell):
    row = int(cell.split(',')[0][1])
    col = int(cell.split(',')[1])
    adjacent = []
    
    # Check bounds and add valid adjacent cells
    if row > 1: adjacent.append(f"C{row-1},{col}")
    if row < 4: adjacent.append(f"C{row+1},{col}")
    if col > 1: adjacent.append(f"C{row},{col-1}")
    if col < 8: adjacent.append(f"C{row},{col+1}")
    
    return adjacent

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def get_next_states(current_state):
    states = []
    # Sort boxes by distance to goal to prioritize moving boxes that are far from their goals
    boxes = sorted(current_state.keys(), 
                  key=lambda b: manhattan_distance(current_state[b], goal_state[b]),
                  reverse=True)
    
    for box in boxes:
        if current_state[box] == goal_state[box]:
            continue
        
        current_pos = current_state[box]
        for next_pos in get_adjacent_cells(current_pos):
            if next_pos not in current_state.values():
                new_state = current_state.copy()
                new_state[box] = next_pos
                states.append(new_state)
    
    return states

def solve():
    visited = set()
    queue = [(0, 0, initial_state, [initial_state])]
    visited.add(str(initial_state))
    
    while queue:
        _, steps, current_state, path = heappop(queue)
        
        if current_state == goal_state:
            return path
        
        for next_state in get_next_states(current_state):
            state_str = str(next_state)
            if state_str not in visited:
                visited.add(state_str)
                h_score = sum(manhattan_distance(next_state[b], goal_state[b]) for b in next_state)
                new_path = path + [next_state]
                heappush(queue, (steps + 1 + h_score, steps + 1, next_state, new_path))
                
                # Early stopping if path gets too long
                if len(new_path) > 50:
                    continue
    
    return None

# Find solution and print
solution = solve()
if solution:
    # Format solution for output
    formatted_solution = []
    for state in solution:
        formatted_solution.append(dict(sorted(state.items())))
    print(json.dumps(formatted_solution))
else:
    print("No solution found")