def generate_plan(initial_state, goal_state, adjacency):
    plan = [initial_state.copy()]
    current_state = initial_state.copy()

    def move_box(box, path):
        for step in path:
            current_state[box] = step
            plan.append(current_state.copy())

    # Define paths for each box
    paths = {
        "box1": ["C1,4", "C1,5", "C1,6", "C1,7", "C1,8", "C2,8"],
        "box2": [],  # Already at goal
        "box3": ["C2,2", "C1,2"],
        "box4": ["C2,1", "C2,2", "C2,3", "C2,4", "C2,5"]
    }

    # Move each box along its path
    for box, path in paths.items():
        move_box(box, path)

    return plan

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

# Generate and print the plan
plan = generate_plan(initial_state, goal_state, adjacency)
print(plan)