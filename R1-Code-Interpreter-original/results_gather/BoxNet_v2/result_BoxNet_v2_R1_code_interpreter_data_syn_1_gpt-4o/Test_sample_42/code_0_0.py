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

# Plan generation
def generate_plan(initial, goal, adjacency):
    plan = [initial.copy()]
    current_state = initial.copy()

    # Function to move a box
    def move_box(box, target):
        nonlocal current_state
        path = [current_state[box]]
        while path[-1] != target:
            for neighbor in adjacency[path[-1]]:
                if neighbor not in current_state.values() or neighbor == target:
                    path.append(neighbor)
                    break
        for step in path[1:]:
            current_state[box] = step
            plan.append(current_state.copy())

    # Move each box to its goal
    for box, target in goal.items():
        if current_state[box] != target:
            move_box(box, target)

    return plan

# Generate the plan
plan = generate_plan(initial_state, goal_state, adjacency)

# Output the plan as JSON
print(json.dumps(plan, indent=2))