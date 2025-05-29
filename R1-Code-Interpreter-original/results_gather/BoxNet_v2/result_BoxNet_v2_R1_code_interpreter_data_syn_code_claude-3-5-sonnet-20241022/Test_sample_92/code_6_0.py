import json

def get_minimal_solution():
    # Start with initial state
    solution = [{
        "box1": "C1,4",
        "box2": "C2,8",
        "box3": "C2,2",
        "box4": "C2,1"
    }]
    
    # Move box1 to goal (C2,8)
    current = dict(solution[-1])
    for pos in ["C1,5", "C1,6", "C1,7", "C1,8", "C2,8"]:
        current = dict(current)
        current["box1"] = pos
        solution.append(current)
    
    # Move box2 to goal (C2,2)
    current = dict(solution[-1])
    for pos in ["C2,7", "C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]:
        current = dict(current)
        current["box2"] = pos
        solution.append(current)
    
    # Move box3 to goal (C1,2)
    current = dict(solution[-1])
    current = dict(current)
    current["box3"] = "C1,2"
    solution.append(current)
    
    # Move box4 to goal (C2,5)
    current = dict(solution[-1])
    for pos in ["C2,2", "C2,3", "C2,4", "C2,5"]:
        current = dict(current)
        current["box4"] = pos
        solution.append(current)
    
    return solution

# Generate solution and print directly with minimal formatting
result = get_minimal_solution()
print(json.dumps(result, separators=(',', ':')))