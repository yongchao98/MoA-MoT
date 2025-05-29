import json

initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
    "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5",
    "box4": "C3,1", "box5": "C1,3"
}

# Hardcoded solution path based on analysis
solution = [
    {"box1": "C5,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},  # Initial state
    {"box1": "C4,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},  # Move box1 toward goal
    {"box1": "C3,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},  # Move box1 toward goal
    {"box1": "C2,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},  # box1 reaches goal
    {"box1": "C2,5", "box2": "C2,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},  # Move box2 toward goal
    {"box1": "C2,5", "box2": "C3,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},  # Move box2 toward goal
    {"box1": "C2,5", "box2": "C3,4", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},  # box2 reaches goal
    {"box1": "C2,5", "box2": "C3,4", "box3": "C2,3", "box4": "C4,5", "box5": "C4,4"},  # Move box3 toward goal
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,3", "box4": "C4,5", "box5": "C4,4"},  # Move box3 toward goal
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,4", "box4": "C4,5", "box5": "C4,4"},  # Move box3 toward goal
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C4,5", "box5": "C4,4"},  # box3 reaches goal
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C4,4", "box5": "C4,4"},  # Move box4 toward goal
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C4,3", "box5": "C4,4"},  # Move box4 toward goal
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,3", "box5": "C4,4"},  # Move box4 toward goal
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,2", "box5": "C4,4"},  # Move box4 toward goal
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C4,4"},  # box4 reaches goal
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C3,4"},  # Move box5 toward goal
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C2,4"},  # Move box5 toward goal
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,4"},  # Move box5 toward goal
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"}   # box5 reaches goal
]

print(json.dumps(solution))