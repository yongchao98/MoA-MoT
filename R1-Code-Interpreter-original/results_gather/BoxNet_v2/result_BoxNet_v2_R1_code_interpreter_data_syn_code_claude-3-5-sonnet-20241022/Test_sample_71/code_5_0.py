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

def solve():
    current_state = initial_state
    path = [current_state]
    
    # Move box4 to C2,5
    while current_state["box4"] != "C2,5":
        if current_state["box4"] == "C2,1":
            next_state = deepcopy(current_state)
            next_state["box4"] = "C2,2"
            path.append(next_state)
            current_state = next_state
        elif current_state["box4"] == "C2,2":
            next_state = deepcopy(current_state)
            next_state["box4"] = "C2,3"
            path.append(next_state)
            current_state = next_state
        elif current_state["box4"] == "C2,3":
            next_state = deepcopy(current_state)
            next_state["box4"] = "C2,4"
            path.append(next_state)
            current_state = next_state
        elif current_state["box4"] == "C2,4":
            next_state = deepcopy(current_state)
            next_state["box4"] = "C2,5"
            path.append(next_state)
            current_state = next_state
    
    # Move box3 to C1,2
    while current_state["box3"] != "C1,2":
        if current_state["box3"] == "C2,2":
            next_state = deepcopy(current_state)
            next_state["box3"] = "C1,2"
            path.append(next_state)
            current_state = next_state
    
    # Move box2 to C2,2
    while current_state["box2"] != "C2,2":
        if current_state["box2"] == "C2,8":
            next_state = deepcopy(current_state)
            next_state["box2"] = "C2,7"
            path.append(next_state)
            current_state = next_state
        elif current_state["box2"] == "C2,7":
            next_state = deepcopy(current_state)
            next_state["box2"] = "C2,6"
            path.append(next_state)
            current_state = next_state
        elif current_state["box2"] == "C2,6":
            next_state = deepcopy(current_state)
            next_state["box2"] = "C2,5"
            path.append(next_state)
            current_state = next_state
        elif current_state["box2"] == "C2,5":
            next_state = deepcopy(current_state)
            next_state["box2"] = "C2,4"
            path.append(next_state)
            current_state = next_state
        elif current_state["box2"] == "C2,4":
            next_state = deepcopy(current_state)
            next_state["box2"] = "C2,3"
            path.append(next_state)
            current_state = next_state
        elif current_state["box2"] == "C2,3":
            next_state = deepcopy(current_state)
            next_state["box2"] = "C2,2"
            path.append(next_state)
            current_state = next_state
    
    # Move box1 to C2,8
    while current_state["box1"] != "C2,8":
        if current_state["box1"] == "C1,4":
            next_state = deepcopy(current_state)
            next_state["box1"] = "C1,5"
            path.append(next_state)
            current_state = next_state
        elif current_state["box1"] == "C1,5":
            next_state = deepcopy(current_state)
            next_state["box1"] = "C1,6"
            path.append(next_state)
            current_state = next_state
        elif current_state["box1"] == "C1,6":
            next_state = deepcopy(current_state)
            next_state["box1"] = "C1,7"
            path.append(next_state)
            current_state = next_state
        elif current_state["box1"] == "C1,7":
            next_state = deepcopy(current_state)
            next_state["box1"] = "C1,8"
            path.append(next_state)
            current_state = next_state
        elif current_state["box1"] == "C1,8":
            next_state = deepcopy(current_state)
            next_state["box1"] = "C2,8"
            path.append(next_state)
            current_state = next_state
    
    return path

# Solve and print result
solution = solve()
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No solution found")