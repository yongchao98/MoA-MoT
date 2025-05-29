import json

initial_state = {
    "box1": "C4,3", "box2": "C4,6", "box3": "C2,4", "box4": "C1,2"
}

goal_state = {
    "box1": "C2,1", "box2": "C1,4", "box3": "C1,6", "box4": "C4,2"
}

# Simplified solution with predetermined moves
solution = [
    {"box1": "C4,3", "box2": "C4,6", "box3": "C2,4", "box4": "C1,2"},  # Initial state
    {"box1": "C4,2", "box2": "C4,6", "box3": "C2,4", "box4": "C1,2"},  # Move box1
    {"box1": "C4,1", "box2": "C4,6", "box3": "C2,4", "box4": "C1,2"},
    {"box1": "C3,1", "box2": "C4,6", "box3": "C2,4", "box4": "C1,2"},
    {"box1": "C2,1", "box2": "C4,6", "box3": "C2,4", "box4": "C1,2"},
    {"box1": "C2,1", "box2": "C4,5", "box3": "C2,4", "box4": "C1,2"},  # Move box2
    {"box1": "C2,1", "box2": "C3,5", "box3": "C2,4", "box4": "C1,2"},
    {"box1": "C2,1", "box2": "C2,5", "box3": "C2,4", "box4": "C1,2"},
    {"box1": "C2,1", "box2": "C1,5", "box3": "C2,4", "box4": "C1,2"},
    {"box1": "C2,1", "box2": "C1,4", "box3": "C2,4", "box4": "C1,2"},
    {"box1": "C2,1", "box2": "C1,4", "box3": "C2,5", "box4": "C1,2"},  # Move box3
    {"box1": "C2,1", "box2": "C1,4", "box3": "C2,6", "box4": "C1,2"},
    {"box1": "C2,1", "box2": "C1,4", "box3": "C1,6", "box4": "C1,2"},
    {"box1": "C2,1", "box2": "C1,4", "box3": "C1,6", "box4": "C2,2"},  # Move box4
    {"box1": "C2,1", "box2": "C1,4", "box3": "C1,6", "box4": "C3,2"},
    {"box1": "C2,1", "box2": "C1,4", "box3": "C1,6", "box4": "C4,2"}   # Final state
]

# Print solution using json.dumps to ensure proper formatting
print(json.dumps(solution))