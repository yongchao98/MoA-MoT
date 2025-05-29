import json
from collections import deque

# Define the initial and goal states
initial_state = {
    "box1": "C4,3",
    "box2": "C4,6",
    "box3": "C2,4",
    "box4": "C1,2"
}

goal_state = {
    "box1": "C2,1",
    "box2": "C1,4",
    "box3": "C1,6",
    "box4": "C4,2"
}

# Define the adjacency list
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"],
    "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"],
    "C4,6": ["C4,5", "C3,6"]
}

# Function to generate the plan
def generate_plan(initial_state, goal_state, adjacency):
    # Initialize the queue with the initial state
    queue = deque([(initial_state, [initial_state])])
    visited = set()
    visited.add(tuple(initial_state.items()))

    while queue:
        current_state, path = queue.popleft()

        # Check if the current state is the goal state
        if current_state == goal_state:
            return path

        # Generate new states by moving each box
        for box, current_cell in current_state.items():
            for adjacent_cell in adjacency[current_cell]:
                # Create a new state by moving the box to the adjacent cell
                new_state = current_state.copy()
                new_state[box] = adjacent_cell

                # Ensure no two boxes occupy the same cell
                if len(set(new_state.values())) == len(new_state):
                    new_state_tuple = tuple(new_state.items())
                    if new_state_tuple not in visited:
                        visited.add(new_state_tuple)
                        queue.append((new_state, path + [new_state]))

# Generate the plan
plan = generate_plan(initial_state, goal_state, adjacency)

# Print the plan in JSON format
print(json.dumps(plan, indent=2))