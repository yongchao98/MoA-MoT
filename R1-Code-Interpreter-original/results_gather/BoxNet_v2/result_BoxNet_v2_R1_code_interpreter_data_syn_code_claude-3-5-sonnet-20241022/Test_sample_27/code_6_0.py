import json

def get_solution():
    states = []
    
    # Initial state
    state = {
        "box1": "C5,5",
        "box2": "C1,3",
        "box3": "C3,3",
        "box4": "C4,5",
        "box5": "C4,4"
    }
    states.append(dict(state))
    
    # Step by step movements
    
    # 1. First move box2 to C3,4
    state["box2"] = "C1,4"; states.append(dict(state))
    state["box2"] = "C2,4"; states.append(dict(state))
    state["box2"] = "C3,4"; states.append(dict(state))
    
    # 2. Move box3 to C1,5
    state["box3"] = "C2,3"; states.append(dict(state))
    state["box3"] = "C2,4"; states.append(dict(state))
    state["box3"] = "C2,5"; states.append(dict(state))
    state["box3"] = "C1,5"; states.append(dict(state))
    
    # 3. Move box5 to C1,3
    state["box5"] = "C4,3"; states.append(dict(state))
    state["box5"] = "C3,3"; states.append(dict(state))
    state["box5"] = "C2,3"; states.append(dict(state))
    state["box5"] = "C1,3"; states.append(dict(state))
    
    # 4. Move box4 to C3,1
    state["box4"] = "C4,4"; states.append(dict(state))
    state["box4"] = "C4,3"; states.append(dict(state))
    state["box4"] = "C4,2"; states.append(dict(state))
    state["box4"] = "C4,1"; states.append(dict(state))
    state["box4"] = "C3,1"; states.append(dict(state))
    
    # 5. Finally move box1 to C2,5
    state["box1"] = "C4,5"; states.append(dict(state))
    state["box1"] = "C3,5"; states.append(dict(state))
    state["box1"] = "C2,5"; states.append(dict(state))
    
    return states

# Generate and print solution
solution = get_solution()
print(json.dumps(solution))