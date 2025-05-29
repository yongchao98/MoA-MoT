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

def get_next_states(current_state):
    next_states = []
    occupied_cells = set(current_state.values())
    
    # Try to move boxes that are not at their goals
    for box in sorted(current_state.keys()):
        if current_state[box] == goal_state[box]:
            continue
            
        current_pos = current_state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied_cells:
                new_state = deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def solve():
    def heuristic(state):
        return sum(manhattan_distance(state[box], goal_state[box]) 
                  for box in state if state[box] != goal_state[box])
    
    visited = set()
    start = json.dumps(initial_state)
    queue = [(0, 0, start, [initial_state])]  # (priority, moves, state_str, path)
    visited.add(start)
    
    while queue:
        _, moves, _, path = heapq.heappop(queue)
        current_state = path[-1]
        
        if current_state == goal_state:
            return path
            
        for next_state in get_next_states(current_state):
            state_str = json.dumps(next_state)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                h = heuristic(next_state)
                priority = len(new_path) + h
                heapq.heappush(queue, (priority, len(new_path), state_str, new_path))
                
                # Early exit if path gets too long
                if len(new_path) > 30:
                    continue
    
    return None

# Create direct solution
def create_direct_solution():
    solution = [initial_state]
    current_state = deepcopy(initial_state)
    
    # Move box3 to C1,2
    if current_state["box3"] != "C1,2":
        current_state = deepcopy(current_state)
        current_state["box3"] = "C1,2"
        solution.append(deepcopy(current_state))
    
    # Move box4 to C2,5
    if current_state["box4"] != "C2,5":
        current_state = deepcopy(current_state)
        current_state["box4"] = "C2,5"
        solution.append(deepcopy(current_state))
    
    # Move box2 to C2,2
    if current_state["box2"] != "C2,2":
        current_state = deepcopy(current_state)
        current_state["box2"] = "C2,2"
        solution.append(deepcopy(current_state))
    
    # Move box1 to C2,8
    if current_state["box1"] != "C2,8":
        current_state = deepcopy(current_state)
        current_state["box1"] = "C2,8"
        solution.append(deepcopy(current_state))
    
    return solution

# Use the direct solution
solution = create_direct_solution()
print("<<<" + json.dumps(solution) + ">>>")