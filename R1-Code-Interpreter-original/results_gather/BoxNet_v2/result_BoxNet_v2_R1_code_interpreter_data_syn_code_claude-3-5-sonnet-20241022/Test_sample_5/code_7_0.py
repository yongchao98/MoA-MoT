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

def solve():
    solution = []
    current = initial_state.copy()
    solution.append(current.copy())
    
    # Step 1: Move box1 to C2,8
    moves_box1 = ["C1,5", "C1,6", "C1,7", "C1,8", "C2,8"]
    for pos in moves_box1:
        current = current.copy()
        current["box1"] = pos
        solution.append(current.copy())
    
    # Step 2: Move box2 to C2,2
    moves_box2 = ["C2,7", "C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]
    for pos in moves_box2:
        current = current.copy()
        current["box2"] = pos
        solution.append(current.copy())
    
    # Step 3: Move box3 to C1,2
    current = current.copy()
    current["box3"] = "C1,2"
    solution.append(current.copy())
    
    # Step 4: Move box4 to C2,5
    moves_box4 = ["C2,2", "C2,3", "C2,4", "C2,5"]
    for pos in moves_box4:
        current = current.copy()
        current["box4"] = pos
        solution.append(current.copy())
    
    return solution

# Generate solution
solution = solve()
print(json.dumps(solution))