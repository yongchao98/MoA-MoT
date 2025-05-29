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

# Simplified solution - move boxes one by one in a predetermined order
def get_path():
    path = [initial_state]
    current = dict(initial_state)
    
    # First move box3 to C1,2
    moves_box3 = ["C1,2"]
    for pos in moves_box3:
        new_state = dict(current)
        new_state["box3"] = pos
        path.append(new_state)
        current = new_state
    
    # Move box4 to C2,5
    moves_box4 = ["C2,2", "C2,3", "C2,4", "C2,5"]
    for pos in moves_box4:
        new_state = dict(current)
        new_state["box4"] = pos
        path.append(new_state)
        current = new_state
    
    # Move box2 to C2,2
    moves_box2 = ["C2,7", "C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]
    for pos in moves_box2:
        new_state = dict(current)
        new_state["box2"] = pos
        path.append(new_state)
        current = new_state
    
    # Finally move box1 to C2,8
    moves_box1 = ["C1,5", "C1,6", "C1,7", "C1,8", "C2,8"]
    for pos in moves_box1:
        new_state = dict(current)
        new_state["box1"] = pos
        path.append(new_state)
        current = new_state
    
    return path

# Generate and print solution
solution = get_path()
print("<<<" + json.dumps(solution) + ">>>")