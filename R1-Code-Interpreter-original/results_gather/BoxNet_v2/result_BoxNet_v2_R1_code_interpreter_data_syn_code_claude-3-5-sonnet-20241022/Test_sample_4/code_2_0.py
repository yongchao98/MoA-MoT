import copy
from collections import deque

initial_state = {
    "box1": "C4,3",
    "box2": "C4,6",
    "box3": "C2,4",
    "box4": "C1,2"
}

goal_state = {
    "box1": "C2,1",
    "box2": "C1,4",
    "box3": "C1,6",
    "box4": "C4,2"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"],
    "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"],
    "C4,6": ["C4,5", "C3,6"]
}

def find_path(start, goal, occupied_cells):
    if start == goal:
        return [start]
    
    visited = set()
    queue = deque([(start, [start])])
    visited.add(start)
    
    while queue:
        current, path = queue.popleft()
        for next_cell in adjacency[current]:
            if next_cell not in visited and next_cell not in occupied_cells:
                if next_cell == goal:
                    return path + [next_cell]
                visited.add(next_cell)
                queue.append((next_cell, path + [next_cell]))
    return None

def move_boxes():
    states = [initial_state]
    current_state = copy.deepcopy(initial_state)
    
    # First move box1 to its goal
    while current_state["box1"] != goal_state["box1"]:
        occupied = {pos for box, pos in current_state.items() if box != "box1"}
        path = find_path(current_state["box1"], goal_state["box1"], occupied)
        if not path or len(path) < 2:
            return None
        current_state = copy.deepcopy(current_state)
        current_state["box1"] = path[1]
        states.append(copy.deepcopy(current_state))
    
    # Then move box2
    while current_state["box2"] != goal_state["box2"]:
        occupied = {pos for box, pos in current_state.items() if box != "box2"}
        path = find_path(current_state["box2"], goal_state["box2"], occupied)
        if not path or len(path) < 2:
            return None
        current_state = copy.deepcopy(current_state)
        current_state["box2"] = path[1]
        states.append(copy.deepcopy(current_state))
    
    # Then move box3
    while current_state["box3"] != goal_state["box3"]:
        occupied = {pos for box, pos in current_state.items() if box != "box3"}
        path = find_path(current_state["box3"], goal_state["box3"], occupied)
        if not path or len(path) < 2:
            return None
        current_state = copy.deepcopy(current_state)
        current_state["box3"] = path[1]
        states.append(copy.deepcopy(current_state))
    
    # Finally move box4
    while current_state["box4"] != goal_state["box4"]:
        occupied = {pos for box, pos in current_state.items() if box != "box4"}
        path = find_path(current_state["box4"], goal_state["box4"], occupied)
        if not path or len(path) < 2:
            return None
        current_state = copy.deepcopy(current_state)
        current_state["box4"] = path[1]
        states.append(copy.deepcopy(current_state))
    
    return states

# Generate and print the plan
plan = move_boxes()
if plan:
    print(str(plan).replace("'", '"'))
else:
    print("No solution found")