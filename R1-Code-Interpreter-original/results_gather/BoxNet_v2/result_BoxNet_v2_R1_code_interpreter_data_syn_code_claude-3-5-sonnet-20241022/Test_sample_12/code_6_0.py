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

def generate_solution():
    # Predefined solution based on analysis of the problem
    solution = []
    current_state = initial_state.copy()
    solution.append(current_state.copy())
    
    # Move box1 to C2,8
    path_box1 = ["C1,4", "C1,5", "C1,6", "C1,7", "C1,8", "C2,8"]
    for pos in path_box1[1:]:  # Skip first position as it's the initial
        current_state = current_state.copy()
        current_state["box1"] = pos
        solution.append(current_state.copy())
    
    # Move box2 to C2,2
    path_box2 = ["C2,8", "C2,7", "C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]
    for pos in path_box2[1:]:
        current_state = current_state.copy()
        current_state["box2"] = pos
        solution.append(current_state.copy())
    
    # Move box3 to C1,2
    path_box3 = ["C2,2", "C1,2"]
    for pos in path_box3[1:]:
        current_state = current_state.copy()
        current_state["box3"] = pos
        solution.append(current_state.copy())
    
    # Move box4 to C2,5
    path_box4 = ["C2,1", "C2,2", "C2,3", "C2,4", "C2,5"]
    for pos in path_box4[1:]:
        current_state = current_state.copy()
        current_state["box4"] = pos
        solution.append(current_state.copy())
    
    return solution

# Generate and print solution
solution = generate_solution()
print(json.dumps(solution))