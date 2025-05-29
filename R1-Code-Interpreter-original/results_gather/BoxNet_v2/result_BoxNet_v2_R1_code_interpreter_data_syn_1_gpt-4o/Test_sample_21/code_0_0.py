import json

# Initial and goal states
initial_state = {
    "box1": "C1,3",
    "box2": "C2,4",
    "box3": "C1,6",
    "box4": "C2,2",
    "box5": "C1,4"
}

goal_state = {
    "box1": "C2,1",
    "box2": "C1,2",
    "box3": "C1,5",
    "box4": "C2,2",
    "box5": "C2,6"
}

# Adjacency list
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"],
    "C2,6": ["C2,5", "C1,6"]
}

# Plan list
plan = [initial_state.copy()]

# Function to find a path for a single box
def move_box(box, current_state, goal_state):
    current_position = current_state[box]
    goal_position = goal_state[box]
    
    if current_position == goal_position:
        return current_state
    
    # Find a valid move
    for neighbor in adjacency[current_position]:
        if neighbor not in current_state.values():
            current_state[box] = neighbor
            return current_state

# Generate the plan
while plan[-1] != goal_state:
    current_state = plan[-1].copy()
    for box in initial_state.keys():
        if current_state[box] != goal_state[box]:
            current_state = move_box(box, current_state, goal_state)
            plan.append(current_state.copy())
            break

# Output the plan
print(json.dumps(plan, indent=2))