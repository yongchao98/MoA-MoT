import json
from copy import deepcopy
from heapq import heappush, heappop

initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"], "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"], "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"], "C2,5": ["C2,4", "C1,5", "C3,5"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"], "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"], "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"], "C4,5": ["C4,4", "C3,5", "C5,5"],
    "C5,1": ["C5,2", "C4,1"], "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"], "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C4,5"]
}

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_next_states(current_state):
    next_states = []
    occupied_cells = set(current_state.values())
    
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
        
        current_pos = current_state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied_cells:
                new_state = current_state.copy()
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def solve():
    states = [initial_state]
    current_state = initial_state.copy()
    
    # First move box2 to its goal
    while current_state["box2"] != "C3,4":
        if current_state["box2"] == "C1,3":
            current_state["box2"] = "C1,4"
        elif current_state["box2"] == "C1,4":
            current_state["box2"] = "C2,4"
        elif current_state["box2"] == "C2,4":
            current_state["box2"] = "C3,4"
        states.append(current_state.copy())
    
    # Move box3 to its goal
    while current_state["box3"] != "C1,5":
        if current_state["box3"] == "C3,3":
            current_state["box3"] = "C2,3"
        elif current_state["box3"] == "C2,3":
            current_state["box3"] = "C2,4"
        elif current_state["box3"] == "C2,4":
            current_state["box3"] = "C2,5"
        elif current_state["box3"] == "C2,5":
            current_state["box3"] = "C1,5"
        states.append(current_state.copy())
    
    # Move box5 to its goal
    while current_state["box5"] != "C1,3":
        if current_state["box5"] == "C4,4":
            current_state["box5"] = "C4,3"
        elif current_state["box5"] == "C4,3":
            current_state["box5"] = "C3,3"
        elif current_state["box5"] == "C3,3":
            current_state["box5"] = "C2,3"
        elif current_state["box5"] == "C2,3":
            current_state["box5"] = "C1,3"
        states.append(current_state.copy())
    
    # Move box4 to its goal
    while current_state["box4"] != "C3,1":
        if current_state["box4"] == "C4,5":
            current_state["box4"] = "C4,4"
        elif current_state["box4"] == "C4,4":
            current_state["box4"] = "C4,3"
        elif current_state["box4"] == "C4,3":
            current_state["box4"] = "C4,2"
        elif current_state["box4"] == "C4,2":
            current_state["box4"] = "C4,1"
        elif current_state["box4"] == "C4,1":
            current_state["box4"] = "C3,1"
        states.append(current_state.copy())
    
    # Finally move box1 to its goal
    while current_state["box1"] != "C2,5":
        if current_state["box1"] == "C5,5":
            current_state["box1"] = "C4,5"
        elif current_state["box1"] == "C4,5":
            current_state["box1"] = "C3,5"
        elif current_state["box1"] == "C3,5":
            current_state["box1"] = "C2,5"
        states.append(current_state.copy())
    
    return states

# Find and print the solution
solution = solve()
print(json.dumps(solution))