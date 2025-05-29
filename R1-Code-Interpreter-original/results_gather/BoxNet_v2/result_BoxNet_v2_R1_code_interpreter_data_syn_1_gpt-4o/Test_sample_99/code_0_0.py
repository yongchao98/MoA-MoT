# Define the initial and goal states
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

# Define the adjacency list
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

# Function to check if a move is valid
def is_valid_move(current_state, box, new_position):
    # Check if the new position is adjacent and not occupied
    current_position = current_state[box]
    if new_position in adjacency[current_position]:
        for other_box, position in current_state.items():
            if other_box != box and position == new_position:
                return False
        return True
    return False

# Plan generation
def generate_plan(initial_state, goal_state):
    current_state = initial_state.copy()
    plan = [current_state.copy()]
    
    while current_state != goal_state:
        for box, goal_position in goal_state.items():
            if current_state[box] != goal_position:
                # Try to move the box towards its goal
                for new_position in adjacency[current_state[box]]:
                    if is_valid_move(current_state, box, new_position):
                        current_state[box] = new_position
                        plan.append(current_state.copy())
                        break
                break  # Move one box at a time
    return plan

# Generate the plan
plan = generate_plan(initial_state, goal_state)
print(plan)