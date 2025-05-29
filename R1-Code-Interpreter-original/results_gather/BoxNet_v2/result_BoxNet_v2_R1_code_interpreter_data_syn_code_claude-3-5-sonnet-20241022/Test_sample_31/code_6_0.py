import json

def get_minimal_solution():
    solution = [
        # Initial state
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"},
        
        # Move box4 to goal (C3,3)
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C3,4", "box5": "C4,7", "box6": "C3,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C3,3", "box5": "C4,7", "box6": "C3,6"},
        
        # Move box5 to goal (C2,3)
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C3,3", "box5": "C3,7", "box6": "C3,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C3,3", "box5": "C2,7", "box6": "C3,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C3,3", "box5": "C2,3", "box6": "C3,6"},
        
        # Move box3 to goal (C2,7)
        {"box1": "C1,7", "box2": "C2,8", "box3": "C3,5", "box4": "C3,3", "box5": "C2,3", "box6": "C3,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C2,5", "box4": "C3,3", "box5": "C2,3", "box6": "C3,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C3,6"},
        
        # Move box2 to goal (C4,1)
        {"box1": "C1,7", "box2": "C3,8", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C3,6"},
        {"box1": "C1,7", "box2": "C4,8", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C3,6"},
        {"box1": "C1,7", "box2": "C4,7", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C3,6"},
        {"box1": "C1,7", "box2": "C4,6", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C4,5", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C4,4", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C4,3", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C4,2", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C4,1", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        
        # Move box1 to goal (C3,6)
        {"box1": "C2,7", "box2": "C4,1", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C3,7", "box2": "C4,1", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C3,6", "box2": "C4,1", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"}
    ]
    return solution

# Get and print the first few states to verify
solution = get_minimal_solution()
first_states = solution[:5]
print(json.dumps(first_states))