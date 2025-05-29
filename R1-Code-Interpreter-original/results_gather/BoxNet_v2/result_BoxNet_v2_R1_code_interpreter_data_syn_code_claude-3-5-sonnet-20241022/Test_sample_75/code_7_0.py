import json

def generate_solution():
    # Define states
    solution = []
    
    # Phase 1: Initial state
    state = {
        "box1": "C1,4",
        "box2": "C2,8",
        "box3": "C2,2",
        "box4": "C2,1"
    }
    solution.append(dict(state))
    
    # Phase 2: Move box3 to C1,2
    state["box3"] = "C1,2"
    solution.append(dict(state))
    
    # Phase 3: Move box4 to C2,5
    for pos in ["C2,2", "C2,3", "C2,4", "C2,5"]:
        state["box4"] = pos
        solution.append(dict(state))
    
    # Phase 4: Move box2 to C2,2
    for pos in ["C2,7", "C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]:
        state["box2"] = pos
        solution.append(dict(state))
    
    # Phase 5: Move box1 to C2,8
    for pos in ["C1,5", "C1,6", "C1,7", "C1,8", "C2,8"]:
        state["box1"] = pos
        solution.append(dict(state))
    
    return solution

# Generate solution and print as JSON
result = generate_solution()
print(json.dumps(result, separators=(',', ':')))