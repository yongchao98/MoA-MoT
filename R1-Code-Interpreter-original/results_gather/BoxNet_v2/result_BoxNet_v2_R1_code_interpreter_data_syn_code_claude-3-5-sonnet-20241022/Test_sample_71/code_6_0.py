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

def solve():
    path = []
    current_state = initial_state.copy()
    path.append(current_state.copy())
    
    # First move box3 to C1,2 since it's blocking box2's goal
    current_state["box3"] = "C1,2"
    path.append(current_state.copy())
    
    # Move box4 to C2,2 temporarily
    current_state["box4"] = "C2,2"
    path.append(current_state.copy())
    
    # Move box4 to C2,3
    current_state["box4"] = "C2,3"
    path.append(current_state.copy())
    
    # Move box4 to C2,4
    current_state["box4"] = "C2,4"
    path.append(current_state.copy())
    
    # Move box4 to final position C2,5
    current_state["box4"] = "C2,5"
    path.append(current_state.copy())
    
    # Move box2 to C2,7
    current_state["box2"] = "C2,7"
    path.append(current_state.copy())
    
    # Move box2 to C2,6
    current_state["box2"] = "C2,6"
    path.append(current_state.copy())
    
    # Move box2 to C2,5
    current_state["box2"] = "C2,4"
    path.append(current_state.copy())
    
    # Move box2 to C2,3
    current_state["box2"] = "C2,3"
    path.append(current_state.copy())
    
    # Move box2 to final position C2,2
    current_state["box2"] = "C2,2"
    path.append(current_state.copy())
    
    # Move box1 to C1,5
    current_state["box1"] = "C1,5"
    path.append(current_state.copy())
    
    # Move box1 to C1,6
    current_state["box1"] = "C1,6"
    path.append(current_state.copy())
    
    # Move box1 to C1,7
    current_state["box1"] = "C1,7"
    path.append(current_state.copy())
    
    # Move box1 to C1,8
    current_state["box1"] = "C1,8"
    path.append(current_state.copy())
    
    # Move box1 to final position C2,8
    current_state["box1"] = "C2,8"
    path.append(current_state.copy())
    
    return path

# Solve and print result
solution = solve()
print(f"<<<{json.dumps(solution)}>>>")