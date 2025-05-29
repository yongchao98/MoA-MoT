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

def create_solution():
    states = []
    current_state = dict(initial_state)
    states.append(dict(current_state))
    
    # First move box1 to C2,8
    moves_box1 = ["C1,5", "C1,6", "C1,7", "C1,8", "C2,8"]
    for move in moves_box1:
        current_state = dict(current_state)
        current_state["box1"] = move
        states.append(dict(current_state))
    
    # Then move box2 to C2,2
    moves_box2 = ["C2,7", "C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]
    for move in moves_box2:
        current_state = dict(current_state)
        current_state["box2"] = move
        states.append(dict(current_state))
    
    # Then move box3 to C1,2
    moves_box3 = ["C1,2"]
    for move in moves_box3:
        current_state = dict(current_state)
        current_state["box3"] = move
        states.append(dict(current_state))
    
    # Finally move box4 to C2,5
    moves_box4 = ["C2,2", "C2,3", "C2,4", "C2,5"]
    for move in moves_box4:
        current_state = dict(current_state)
        current_state["box4"] = move
        states.append(dict(current_state))
    
    return states

# Generate and print solution
solution = create_solution()
print(json.dumps(solution))