import json

# Initial and goal states
initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}
goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

def create_path():
    solution = []
    current_state = dict(initial_state)
    solution.append(dict(current_state))
    
    # Step 1: Move box1 to C3,6 (goal)
    moves_box1 = ["C2,7", "C3,7", "C3,6"]
    for move in moves_box1:
        current_state = dict(solution[-1])
        current_state["box1"] = move
        solution.append(current_state)
    
    # Step 2: Move box6 to C2,6 (goal)
    current_state = dict(solution[-1])
    current_state["box6"] = "C2,6"
    solution.append(current_state)
    
    # Step 3: Move box4 to C3,3 (goal)
    moves_box4 = ["C3,4", "C3,3"]
    for move in moves_box4:
        current_state = dict(solution[-1])
        current_state["box4"] = move
        solution.append(current_state)
    
    # Step 4: Move box5 to C2,3 (goal)
    moves_box5 = ["C4,6", "C4,5", "C4,4", "C4,3", "C3,3", "C2,3"]
    for move in moves_box5:
        current_state = dict(solution[-1])
        current_state["box5"] = move
        solution.append(current_state)
    
    # Step 5: Move box3 to C2,7 (goal)
    moves_box3 = ["C3,5", "C2,5", "C2,6", "C2,7"]
    for move in moves_box3:
        current_state = dict(solution[-1])
        current_state["box3"] = move
        solution.append(current_state)
    
    # Step 6: Move box2 to C4,1 (goal)
    moves_box2 = ["C3,8", "C4,8", "C4,7", "C4,6", "C4,5", "C4,4", "C4,3", "C4,2", "C4,1"]
    for move in moves_box2:
        current_state = dict(solution[-1])
        current_state["box2"] = move
        solution.append(current_state)
    
    return solution

# Generate and print solution
solution = create_path()
print(json.dumps(solution))