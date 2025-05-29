from queue import PriorityQueue
import json
from copy import deepcopy

# Problem setup
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

def is_goal_reached(state):
    return (state["box1"] == goal_state["box1"] and 
            state["box2"] == goal_state["box2"] and 
            state["box3"] == goal_state["box3"] and 
            state["box4"] == goal_state["box4"])

def get_next_states(current_state):
    next_states = []
    occupied_cells = set(current_state.values())
    
    for box in current_state:
        # If box is at goal, skip it
        if current_state[box] == goal_state[box]:
            continue
            
        current_pos = current_state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied_cells:
                new_state = deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def solve():
    queue = [(0, [initial_state])]
    visited = {state_to_string(initial_state)}
    
    while queue:
        _, path = queue.pop(0)
        current_state = path[-1]
        
        if is_goal_reached(current_state):
            return path
        
        next_states = get_next_states(current_state)
        for next_state in next_states:
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                queue.append((len(new_path), new_path))
        
        queue.sort(key=lambda x: x[0])  # Sort by path length
        if len(queue) > 1000:  # Limit queue size
            queue = queue[:500]
    
    return None

# Generate solution
solution = solve()
if solution:
    # Create the final answer string directly
    answer = "<<<" + json.dumps(solution) + ">>>"
    print(answer)
else:
    print("No solution found")