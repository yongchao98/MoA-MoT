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
def generate_plan(initial_state, goal_state, adjacency):
    plan = [initial_state.copy()]
    current_state = initial_state.copy()

    # Function to move a box
    def move_box(box, current_pos, goal_pos):
        if current_pos == goal_pos:
            return current_pos
        for neighbor in adjacency[current_pos]:
            if neighbor == goal_pos or neighbor not in current_state.values():
                return neighbor
        return current_pos

    # Generate the plan
    while current_state != goal_state:
        new_state = current_state.copy()
        for box, goal_pos in goal_state.items():
            current_pos = current_state[box]
            if current_pos != goal_pos:
                new_state[box] = move_box(box, current_pos, goal_pos)
        plan.append(new_state)
        current_state = new_state

    return plan

# Generate and print the plan
plan = generate_plan(initial_state, goal_state, adjacency)
print(json.dumps(plan, indent=2))