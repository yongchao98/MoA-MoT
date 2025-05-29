import json

# Create minimal solution with 10 states
solution = [
    # Initial state
    {"box1": "C5,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},
    
    # Move box1 to goal (C2,5)
    {"box1": "C4,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},
    {"box1": "C3,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},
    {"box1": "C2,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},
    
    # Move box4 and box5
    {"box1": "C2,5", "box2": "C1,3", "box3": "C3,3", "box4": "C3,5", "box5": "C3,4"},
    {"box1": "C2,5", "box2": "C1,3", "box3": "C3,3", "box4": "C3,4", "box5": "C3,3"},
    {"box1": "C2,5", "box2": "C1,3", "box3": "C2,3", "box4": "C3,3", "box5": "C3,2"},
    {"box1": "C2,5", "box2": "C1,3", "box3": "C2,3", "box4": "C3,1", "box5": "C2,2"},
    
    # Move box2 and box3
    {"box1": "C2,5", "box2": "C2,3", "box3": "C1,3", "box4": "C3,1", "box5": "C2,2"},
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"}
]

print(json.dumps(solution))