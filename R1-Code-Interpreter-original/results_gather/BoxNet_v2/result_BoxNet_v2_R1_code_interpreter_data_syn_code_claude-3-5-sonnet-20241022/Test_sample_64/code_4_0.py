import json
from copy import deepcopy

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
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"], "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"], "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"], "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"], "C2,8": ["C2,7", "C1,8"]
}

def get_next_states(state):
    next_states = []
    occupied = set(state.values())
    
    for box in state:
        # Skip if box is at goal
        if state[box] == goal_state[box]:
            continue
            
        curr_pos = state[box]
        for next_pos in adjacency[curr_pos]:
            if next_pos not in occupied:
                new_state = deepcopy(state)
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def solve():
    solution = [initial_state]
    current = initial_state
    visited = {json.dumps(initial_state, sort_keys=True)}
    
    # Move box1 to C2,8
    while current["box1"] != "C2,8":
        for next_state in get_next_states(current):
            if json.dumps(next_state, sort_keys=True) not in visited:
                if next_state["box1"] in ["C1,5", "C1,6", "C1,7", "C1,8", "C2,8"]:
                    current = next_state
                    solution.append(current)
                    visited.add(json.dumps(current, sort_keys=True))
                    break
    
    # Move box2 to C2,2
    while current["box2"] != "C2,2":
        for next_state in get_next_states(current):
            if json.dumps(next_state, sort_keys=True) not in visited:
                if next_state["box2"] in ["C2,7", "C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]:
                    current = next_state
                    solution.append(current)
                    visited.add(json.dumps(current, sort_keys=True))
                    break
    
    # Move box3 to C1,2
    while current["box3"] != "C1,2":
        for next_state in get_next_states(current):
            if json.dumps(next_state, sort_keys=True) not in visited:
                if next_state["box3"] == "C1,2":
                    current = next_state
                    solution.append(current)
                    visited.add(json.dumps(current, sort_keys=True))
                    break
    
    # Move box4 to C2,5
    while current["box4"] != "C2,5":
        for next_state in get_next_states(current):
            if json.dumps(next_state, sort_keys=True) not in visited:
                if next_state["box4"] in ["C2,2", "C2,3", "C2,4", "C2,5"]:
                    current = next_state
                    solution.append(current)
                    visited.add(json.dumps(current, sort_keys=True))
                    break
    
    return solution

# Solve and print result
solution = solve()
print(json.dumps(solution))