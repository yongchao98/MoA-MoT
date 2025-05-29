# Define the initial state, goal state, and adjacency list
initial_state = {
    "box1": "C5,5",
    "box2": "C1,3",
    "box3": "C3,3",
    "box4": "C4,5",
    "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5",
    "box2": "C3,4",
    "box3": "C1,5",
    "box4": "C3,1",
    "box5": "C1,3"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C1,5", "C3,5"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"],
    "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C3,5", "C5,5"],
    "C5,1": ["C5,2", "C4,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C4,5"]
}

# Initialize the plan with the initial state
plan = [initial_state.copy()]

# Function to find the next move for a box
def find_next_move(current, goal, occupied):
    if current == goal:
        return current
    for neighbor in adjacency[current]:
        if neighbor not in occupied:
            return neighbor
    return current

# Generate the plan
while plan[-1] != goal_state:
    current_state = plan[-1].copy()
    occupied_cells = set(current_state.values())
    for box, current_location in current_state.items():
        if current_location != goal_state[box]:
            next_location = find_next_move(current_location, goal_state[box], occupied_cells)
            current_state[box] = next_location
            occupied_cells.add(next_location)
    plan.append(current_state)

# Print the plan
print(plan)