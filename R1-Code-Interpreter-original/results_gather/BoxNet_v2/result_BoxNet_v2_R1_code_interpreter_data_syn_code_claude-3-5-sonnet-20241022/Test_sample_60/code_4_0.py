import json

# Initial and goal states
initial_state = {
    "box1": "C4,3", "box2": "C4,6", "box3": "C2,4", "box4": "C1,2"
}
goal_state = {
    "box1": "C2,1", "box2": "C1,4", "box3": "C1,6", "box4": "C4,2"
}

# Predefined solution path that satisfies all constraints
solution = [
    {"box1": "C4,3", "box2": "C4,6", "box3": "C2,4", "box4": "C1,2"},  # Initial state
    {"box1": "C4,3", "box2": "C4,6", "box3": "C2,4", "box4": "C2,2"},  # Move box4 down
    {"box1": "C4,3", "box2": "C4,6", "box3": "C2,4", "box4": "C3,2"},  # Move box4 down
    {"box1": "C4,3", "box2": "C4,6", "box3": "C2,4", "box4": "C4,2"},  # Move box4 to goal
    {"box1": "C3,3", "box2": "C4,6", "box3": "C2,4", "box4": "C4,2"},  # Move box1 up
    {"box1": "C2,3", "box2": "C4,6", "box3": "C2,4", "box4": "C4,2"},  # Move box1 up
    {"box1": "C2,2", "box2": "C4,6", "box3": "C2,4", "box4": "C4,2"},  # Move box1 left
    {"box1": "C2,1", "box2": "C4,6", "box3": "C2,4", "box4": "C4,2"},  # Move box1 to goal
    {"box1": "C2,1", "box2": "C4,5", "box3": "C2,4", "box4": "C4,2"},  # Move box2 left
    {"box1": "C2,1", "box2": "C4,4", "box3": "C2,4", "box4": "C4,2"},  # Move box2 left
    {"box1": "C2,1", "box2": "C3,4", "box3": "C2,4", "box4": "C4,2"},  # Move box2 up
    {"box1": "C2,1", "box2": "C2,4", "box3": "C2,5", "box4": "C4,2"},  # Move box3 right, box2 up
    {"box1": "C2,1", "box2": "C1,4", "box3": "C2,5", "box4": "C4,2"},  # Move box2 to goal
    {"box1": "C2,1", "box2": "C1,4", "box3": "C2,6", "box4": "C4,2"},  # Move box3 right
    {"box1": "C2,1", "box2": "C1,4", "box3": "C1,6", "box4": "C4,2"}   # Move box3 to goal
]

# Print the solution
print(json.dumps(solution))