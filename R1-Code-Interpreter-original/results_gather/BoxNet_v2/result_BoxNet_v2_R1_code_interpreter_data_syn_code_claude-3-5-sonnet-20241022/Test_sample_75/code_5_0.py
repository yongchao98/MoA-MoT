from queue import PriorityQueue
import json
import copy

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

def move_box(current_state, box, path):
    states = []
    state = copy.deepcopy(current_state)
    
    for pos in path:
        new_state = copy.deepcopy(state)
        new_state[box] = pos
        states.append(new_state)
        state = new_state
    
    return states

def find_solution():
    solution = [initial_state]
    current_state = copy.deepcopy(initial_state)
    
    # First move box3 to C1,2
    if current_state["box3"] != "C1,2":
        solution.extend(move_box(current_state, "box3", ["C1,2"]))
        current_state = solution[-1]
    
    # Then move box4 to C2,5
    if current_state["box4"] != "C2,5":
        path = ["C2,2", "C2,3", "C2,4", "C2,5"]
        solution.extend(move_box(current_state, "box4", path))
        current_state = solution[-1]
    
    # Move box2 to C2,2
    if current_state["box2"] != "C2,2":
        path = ["C2,7", "C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]
        solution.extend(move_box(current_state, "box2", path))
        current_state = solution[-1]
    
    # Finally move box1 to C2,8
    if current_state["box1"] != "C2,8":
        path = ["C1,5", "C1,6", "C1,7", "C1,8", "C2,8"]
        solution.extend(move_box(current_state, "box1", path))
    
    return solution

# Generate and print solution
solution = find_solution()
print(json.dumps(solution))