import heapq
from copy import deepcopy
import json

initial_state = {
    "box1": "C1,4",
    "box2": "C2,8",
    "box3": "C2,2",
    "box4": "C2,1"
}

goal_state = {
    "box1": "C2,8",
    "box2": "C2,2",
    "box3": "C1,2",
    "box4": "C2,5"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"],
    "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"],
    "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"],
    "C2,8": ["C2,7", "C1,8"]
}

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def get_next_states(current_state):
    next_states = []
    occupied_cells = set(current_state.values())
    
    # First try to move boxes that are blocking others
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
        
        # Prioritize moving boxes that are in the way of others
        current_pos = current_state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied_cells:
                new_state = deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def solve():
    def heuristic(state):
        total = 0
        for box in state:
            if state[box] != goal_state[box]:
                total += manhattan_distance(state[box], goal_state[box])
        return total

    visited = set()
    queue = [(0, 0, initial_state, [initial_state])]  # (priority, moves, state, path)
    visited.add(state_to_string(initial_state))
    
    while queue:
        _, moves, current_state, path = heapq.heappop(queue)
        
        if current_state == goal_state:
            return path
        
        for next_state in get_next_states(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                h = heuristic(next_state)
                priority = moves + 1 + h
                heapq.heappush(queue, (priority, moves + 1, next_state, new_path))
    
    return None

# Solve and print result
solution = solve()
if solution:
    # Format solution to be more compact
    compact_solution = []
    prev_state = None
    for state in solution:
        if state != prev_state:
            compact_solution.append(state)
            prev_state = state
    print("<<<" + json.dumps(compact_solution) + ">>>")
else:
    print("No solution found")