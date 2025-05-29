import json

initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
    "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5",
    "box4": "C3,1", "box5": "C1,3"
}

# Pre-calculated sequence of moves
solution = [
    {"box1": "C5,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},  # Initial state
    {"box1": "C5,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C3,4"},  # Move box5 up
    {"box1": "C5,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C2,4"},  # Move box5 up
    {"box1": "C5,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C1,4"},  # Move box5 up
    {"box1": "C5,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C1,3"},  # Move box5 left to goal
    {"box1": "C5,5", "box2": "C2,3", "box3": "C3,3", "box4": "C4,5", "box5": "C1,3"},  # Move box2 up
    {"box1": "C5,5", "box2": "C3,3", "box3": "C2,3", "box4": "C4,5", "box5": "C1,3"},  # Move box2 up, box3 down
    {"box1": "C5,5", "box2": "C3,4", "box3": "C2,3", "box4": "C4,5", "box5": "C1,3"},  # Move box2 right to goal
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,3", "box4": "C4,5", "box5": "C1,4"},  # Move box3 up, box5 right
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,4", "box4": "C4,5", "box5": "C1,3"},  # Move box3 right, box5 left
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C4,5", "box5": "C1,3"},  # Move box3 right to goal
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,5", "box5": "C1,3"},  # Move box4 left
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,4", "box5": "C1,3"},  # Move box4 down
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,3", "box5": "C1,3"},  # Move box4 down
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,2", "box5": "C1,3"},  # Move box4 down
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"},  # Move box4 down to goal
    {"box1": "C4,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"},  # Move box1 up
    {"box1": "C3,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"},  # Move box1 up
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"}   # Move box1 up to goal
]

# Verify solution
def verify_solution(solution):
    # Check if first state matches initial state
    if solution[0] != initial_state:
        return False
    
    # Check if last state matches goal state
    if solution[-1] != goal_state:
        return False
    
    # Check each transition
    for i in range(len(solution) - 1):
        current = solution[i]
        next_state = solution[i + 1]
        
        # Count differences
        differences = sum(1 for box in current if current[box] != next_state[box])
        if differences != 1:
            return False
            
    return True

# Verify and print solution
if verify_solution(solution):
    print(json.dumps(solution))
else:
    print("Invalid solution")