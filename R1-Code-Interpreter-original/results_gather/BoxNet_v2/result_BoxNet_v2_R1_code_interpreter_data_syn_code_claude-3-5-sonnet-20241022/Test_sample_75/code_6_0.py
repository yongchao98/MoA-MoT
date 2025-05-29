import json

def solve_boxes():
    # Initial and goal states
    initial_state = {
        "box1": "C1,4", "box2": "C2,8",
        "box3": "C2,2", "box4": "C2,1"
    }
    goal_state = {
        "box1": "C2,8", "box2": "C2,2",
        "box3": "C1,2", "box4": "C2,5"
    }
    
    # Build solution step by step
    solution = []
    current = dict(initial_state)
    solution.append(current.copy())
    
    # Move box3 to C1,2
    current["box3"] = "C1,2"
    solution.append(current.copy())
    
    # Move box4 to goal through C2,2->C2,3->C2,4->C2,5
    moves4 = ["C2,2", "C2,3", "C2,4", "C2,5"]
    for move in moves4:
        current["box4"] = move
        solution.append(current.copy())
    
    # Move box2 from C2,8 to C2,2
    moves2 = ["C2,7", "C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]
    for move in moves2:
        current["box2"] = move
        solution.append(current.copy())
    
    # Move box1 to C2,8
    moves1 = ["C1,5", "C1,6", "C1,7", "C1,8", "C2,8"]
    for move in moves1:
        current["box1"] = move
        solution.append(current.copy())
    
    return solution

# Generate and print solution
solution = solve_boxes()
print(json.dumps(solution))