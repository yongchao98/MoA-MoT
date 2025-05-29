import json

def create_final_solution():
    states = []
    
    # Initial state
    state = {
        "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
        "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
    }
    states.append(state.copy())
    
    # Move box1 to goal (C3,6) and box6 to goal (C2,6)
    state = state.copy()
    state["box1"] = "C2,7"
    states.append(state.copy())
    
    state = state.copy()
    state["box1"] = "C3,7"
    states.append(state.copy())
    
    state = state.copy()
    state["box1"] = "C3,6"
    states.append(state.copy())
    
    state = state.copy()
    state["box6"] = "C2,6"
    states.append(state.copy())
    
    # Move box3 to goal (C2,7)
    state = state.copy()
    state["box3"] = "C3,5"
    states.append(state.copy())
    
    state = state.copy()
    state["box3"] = "C2,5"
    states.append(state.copy())
    
    state = state.copy()
    state["box3"] = "C2,7"
    states.append(state.copy())
    
    # Move box5 to goal (C2,3)
    state = state.copy()
    state["box5"] = "C3,7"
    states.append(state.copy())
    
    state = state.copy()
    state["box5"] = "C2,7"
    states.append(state.copy())
    
    return states

solution = create_final_solution()
print(json.dumps(solution))