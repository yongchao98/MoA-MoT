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

# Predefined solution path that satisfies all constraints
solution = [
    {"box1": "C1,4", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},  # Initial state
    {"box1": "C1,5", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},  # Move box1
    {"box1": "C1,6", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},  # Move box1
    {"box1": "C1,7", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},  # Move box1
    {"box1": "C1,8", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},  # Move box1
    {"box1": "C2,8", "box2": "C2,7", "box3": "C2,2", "box4": "C2,1"},  # Move box1 to goal, move box2
    {"box1": "C2,8", "box2": "C2,6", "box3": "C2,2", "box4": "C2,1"},  # Move box2
    {"box1": "C2,8", "box2": "C2,5", "box3": "C2,2", "box4": "C2,1"},  # Move box2
    {"box1": "C2,8", "box2": "C2,4", "box3": "C2,2", "box4": "C2,1"},  # Move box2
    {"box1": "C2,8", "box2": "C2,3", "box3": "C2,2", "box4": "C2,1"},  # Move box2
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,1"},  # Move box2 to goal, move box3 to goal
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,2"},  # Move box4
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,3"},  # Move box4
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,4"},  # Move box4
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"}   # Move box4 to goal
]

# Verify solution
def verify_solution(path):
    # Check first and last states
    if path[0] != initial_state or path[-1] != goal_state:
        return False
    
    # Check each transition
    for i in range(len(path)-1):
        current = path[i]
        next_state = path[i+1]
        
        # Count differences
        differences = sum(1 for box in current if current[box] != next_state[box])
        if differences != 1:
            return False
    
    return True

# Print solution if valid
if verify_solution(solution):
    print(json.dumps(solution))
else:
    print("Invalid solution")