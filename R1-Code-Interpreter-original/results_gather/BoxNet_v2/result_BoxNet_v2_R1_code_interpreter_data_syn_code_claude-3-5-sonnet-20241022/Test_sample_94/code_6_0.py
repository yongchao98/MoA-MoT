import json

# Create a minimal solution path focusing on essential moves
solution = [
    # Initial state
    {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"},
    # Move box4 to goal (C2,4 -> C3,3)
    {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C3,3", "box5": "C4,7", "box6": "C3,6"},
    # Move box5 to goal (C4,7 -> C2,3)
    {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C3,3", "box5": "C3,7", "box6": "C3,6"},
    {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C3,3", "box5": "C2,7", "box6": "C3,6"},
    {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C3,3", "box5": "C2,3", "box6": "C3,6"},
    # Move box3 to goal (C4,5 -> C2,7)
    {"box1": "C1,7", "box2": "C2,8", "box3": "C3,5", "box4": "C3,3", "box5": "C2,3", "box6": "C3,6"},
    {"box1": "C1,7", "box2": "C2,8", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C3,6"},
    # Move box6 to goal (C3,6 -> C2,6)
    {"box1": "C1,7", "box2": "C2,8", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
    # Move box1 to goal (C1,7 -> C3,6)
    {"box1": "C2,7", "box2": "C2,8", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
    {"box1": "C3,6", "box2": "C2,8", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
    # Move box2 to goal (C2,8 -> C4,1)
    {"box1": "C3,6", "box2": "C3,8", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
    {"box1": "C3,6", "box2": "C4,8", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
    {"box1": "C3,6", "box2": "C4,1", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"}
]

# Convert to string with minimal formatting
result = json.dumps(solution, separators=(',',':'))[:1000]
print(result)