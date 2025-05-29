import json

# Initial and goal states
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

# Adjacency list
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"],
    "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"],
    "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"],
    "C2,8": ["C2,7", "C1,8"]
}

# Plan generation
plan = [initial_state.copy()]

# Function to move a box
def move_box(box, current_state, target):
    if current_state[box] == target:
        return current_state
    for adj in adjacency[current_state[box]]:
        if adj == target or (adj not in current_state.values()):
            new_state = current_state.copy()
            new_state[box] = adj
            return new_state
    return current_state

# Generate the plan
current_state = initial_state.copy()
while current_state != goal_state:
    for box in initial_state.keys():
        if current_state[box] != goal_state[box]:
            new_state = move_box(box, current_state, goal_state[box])
            if new_state != current_state:
                plan.append(new_state)
                current_state = new_state
                break

# Output the plan
print(json.dumps(plan, indent=2))