import json

solution = [
    # Initial state
    {"box1": "C5,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},
    # Move box2 to C3,4
    {"box1": "C5,5", "box2": "C2,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},
    {"box1": "C5,5", "box2": "C3,3", "box3": "C2,3", "box4": "C4,5", "box5": "C4,4"},
    {"box1": "C5,5", "box2": "C3,4", "box3": "C2,3", "box4": "C4,5", "box5": "C4,4"},
    # Move box3 to C1,5
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,3", "box4": "C4,5", "box5": "C4,4"},
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,4", "box4": "C4,5", "box5": "C4,4"},
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C4,5", "box5": "C4,4"},
    # Move box4 to C3,1
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C4,4", "box5": "C4,3"},
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,4", "box5": "C4,3"},
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C4,3"},
    # Move box5 to C1,3
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"},
    # Move box1 to C2,5
    {"box1": "C3,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"},
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"}
]

print(json.dumps(solution))