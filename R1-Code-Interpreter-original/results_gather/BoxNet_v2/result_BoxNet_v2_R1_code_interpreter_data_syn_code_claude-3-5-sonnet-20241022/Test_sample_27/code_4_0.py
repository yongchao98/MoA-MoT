import json

initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"
}

def generate_solution():
    solution = []
    current_state = initial_state.copy()
    solution.append(current_state.copy())
    
    # Phase 1: Move box2 to C3,4
    moves_box2 = ["C1,4", "C2,4", "C3,4"]
    for move in moves_box2:
        current_state = current_state.copy()
        current_state["box2"] = move
        solution.append(current_state.copy())
    
    # Phase 2: Move box3 to C1,5
    moves_box3 = ["C2,3", "C2,4", "C2,5", "C1,5"]
    for move in moves_box3:
        current_state = current_state.copy()
        current_state["box3"] = move
        solution.append(current_state.copy())
    
    # Phase 3: Move box5 to C1,3
    moves_box5 = ["C4,3", "C3,3", "C2,3", "C1,3"]
    for move in moves_box5:
        current_state = current_state.copy()
        current_state["box5"] = move
        solution.append(current_state.copy())
    
    # Phase 4: Move box4 to C3,1
    moves_box4 = ["C4,4", "C4,3", "C4,2", "C4,1", "C3,1"]
    for move in moves_box4:
        current_state = current_state.copy()
        current_state["box4"] = move
        solution.append(current_state.copy())
    
    # Phase 5: Move box1 to C2,5
    moves_box1 = ["C4,5", "C3,5", "C2,5"]
    for move in moves_box1:
        current_state = current_state.copy()
        current_state["box1"] = move
        solution.append(current_state.copy())
    
    return solution

# Generate and print the solution
solution = generate_solution()
print(json.dumps(solution))