from collections import deque

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

# Function to check if a state is valid
def is_valid_state(state):
    # Ensure no two boxes are in the same cell
    locations = list(state.values())
    return len(locations) == len(set(locations))

# Function to generate the plan
def generate_plan(initial_state, goal_state, adjacency):
    # Initialize the queue with the initial state
    queue = deque([(initial_state, [initial_state])])
    visited = set()
    visited.add(tuple(initial_state.items()))

    while queue:
        current_state, path = queue.popleft()

        # Check if we reached the goal state
        if current_state == goal_state:
            return path

        # Try moving each box to its adjacent cells
        for box, current_location in current_state.items():
            for adjacent in adjacency[current_location]:
                new_state = current_state.copy()
                new_state[box] = adjacent

                # Check if the new state is valid and not visited
                if is_valid_state(new_state) and tuple(new_state.items()) not in visited:
                    visited.add(tuple(new_state.items()))
                    queue.append((new_state, path + [new_state]))

# Generate the plan
plan = generate_plan(initial_state, goal_state, adjacency)
print(plan)