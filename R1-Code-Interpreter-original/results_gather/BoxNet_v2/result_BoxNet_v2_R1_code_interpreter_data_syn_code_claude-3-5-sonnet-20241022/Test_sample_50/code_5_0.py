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

def create_solution():
    solution = []
    current_state = copy.deepcopy(initial_state)
    solution.append(current_state)
    
    # First move box1 to C2,8
    moves_box1 = ["C1,5", "C1,6", "C1,7", "C1,8", "C2,8"]
    for move in moves_box1:
        new_state = copy.deepcopy(current_state)
        new_state["box1"] = move
        solution.append(new_state)
        current_state = copy.deepcopy(new_state)
    
    # Then move box2 to C2,2
    moves_box2 = ["C2,7", "C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]
    for move in moves_box2:
        new_state = copy.deepcopy(current_state)
        new_state["box2"] = move
        solution.append(new_state)
        current_state = copy.deepcopy(new_state)
    
    # Then move box3 to C1,2
    moves_box3 = ["C1,2"]
    for move in moves_box3:
        new_state = copy.deepcopy(current_state)
        new_state["box3"] = move
        solution.append(new_state)
        current_state = copy.deepcopy(new_state)
    
    # Finally move box4 to C2,5
    moves_box4 = ["C2,2", "C2,3", "C2,4", "C2,5"]
    for move in moves_box4:
        new_state = copy.deepcopy(current_state)
        new_state["box4"] = move
        solution.append(new_state)
        current_state = copy.deepcopy(new_state)
    
    return solution

# Generate and print solution
solution = create_solution()
print(json.dumps(solution))