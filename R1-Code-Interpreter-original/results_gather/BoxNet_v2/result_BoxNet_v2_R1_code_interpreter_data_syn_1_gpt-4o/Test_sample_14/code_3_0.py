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

# Function to generate the plan
def generate_plan(initial_state, goal_state, adjacency):
    plan = [initial_state.copy()]
    current_state = initial_state.copy()

    while current_state != goal_state:
        moved = False
        for box, current_pos in current_state.items():
            if current_pos != goal_state[box]:
                # Move the box towards its goal
                for neighbor in adjacency[current_pos]:
                    if neighbor == goal_state[box] or (neighbor not in current_state.values()):
                        current_state[box] = neighbor
                        moved = True
                        break
            if moved:
                break
        plan.append(current_state.copy())

    return plan

# Generate the plan
plan = generate_plan(initial_state, goal_state, adjacency)

# Print the plan in JSON format
print(json.dumps(plan, indent=2))