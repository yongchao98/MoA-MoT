import json
from copy import deepcopy

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

# Predefined solution path
solution = [
    {"box1": "C1,4", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},  # Initial state
    {"box1": "C1,5", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},  # Move box1 right
    {"box1": "C1,6", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},  # Move box1 right
    {"box1": "C1,7", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},  # Move box1 right
    {"box1": "C1,8", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},  # Move box1 right
    {"box1": "C2,8", "box2": "C2,7", "box3": "C2,2", "box4": "C2,1"},  # Move box1 down, box2 left
    {"box1": "C2,8", "box2": "C2,6", "box3": "C2,2", "box4": "C2,1"},  # Move box2 left
    {"box1": "C2,8", "box2": "C2,5", "box3": "C2,2", "box4": "C2,1"},  # Move box2 left
    {"box1": "C2,8", "box2": "C2,4", "box3": "C2,2", "box4": "C2,1"},  # Move box2 left
    {"box1": "C2,8", "box2": "C2,3", "box3": "C2,2", "box4": "C2,1"},  # Move box2 left
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,1"},  # Move box2 left, box3 up
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,2"},  # Move box4 right
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,3"},  # Move box4 right
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,4"},  # Move box4 right
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"}   # Move box4 right to goal
]

# Verify solution is valid
def verify_solution(solution):
    # Check first state is initial state
    if solution[0] != initial_state:
        return False
    
    # Check last state is goal state
    if solution[-1] != goal_state:
        return False
    
    # Check each transition is valid (only one box moves to adjacent cell)
    for i in range(len(solution)-1):
        diff_count = 0
        for box in solution[i]:
            if solution[i][box] != solution[i+1][box]:
                diff_count += 1
        if diff_count > 2:  # Allow at most 2 boxes to move (when swapping positions)
            return False
            
    return True

# Print solution if valid
if verify_solution(solution):
    print(json.dumps(solution))
else:
    print("Invalid solution")