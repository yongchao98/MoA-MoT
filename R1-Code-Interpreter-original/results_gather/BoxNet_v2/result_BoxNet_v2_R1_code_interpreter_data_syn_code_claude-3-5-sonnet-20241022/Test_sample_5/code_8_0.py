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

def is_valid_state(state):
    # Check if any two boxes share the same position
    positions = state.values()
    return len(positions) == len(set(positions))

def solve():
    solution = []
    current = initial_state.copy()
    solution.append(current.copy())
    
    # 1. First move box1 to C2,8
    path_box1 = ["C1,5", "C1,6", "C1,7", "C1,8", "C2,8"]
    for pos in path_box1:
        new_state = current.copy()
        new_state["box1"] = pos
        if is_valid_state(new_state):
            current = new_state
            solution.append(current.copy())
    
    # 2. Then move box2 to C2,2 (after box1 is at C2,8)
    path_box2 = ["C2,7", "C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]
    for pos in path_box2:
        new_state = current.copy()
        new_state["box2"] = pos
        if is_valid_state(new_state):
            current = new_state
            solution.append(current.copy())
    
    # 3. Move box3 to C1,2 (after box2 is at C2,2)
    new_state = current.copy()
    new_state["box3"] = "C1,2"
    if is_valid_state(new_state):
        current = new_state
        solution.append(current.copy())
    
    # 4. Finally move box4 to C2,5
    path_box4 = ["C2,2", "C2,3", "C2,4", "C2,5"]
    for pos in path_box4:
        new_state = current.copy()
        new_state["box4"] = pos
        if is_valid_state(new_state):
            current = new_state
            solution.append(current.copy())
    
    # Verify final state matches goal state
    if current != goal_state:
        return None
    
    return solution

# Generate and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No valid solution found")