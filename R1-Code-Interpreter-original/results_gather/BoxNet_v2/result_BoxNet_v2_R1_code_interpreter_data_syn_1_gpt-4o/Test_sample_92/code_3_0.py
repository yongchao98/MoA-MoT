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

# Function to find a path from start to goal using BFS
def find_path(start, goal, adjacency):
    from collections import deque
    queue = deque([(start, [start])])
    visited = set()
    while queue:
        current, path = queue.popleft()
        if current == goal:
            return path
        visited.add(current)
        for neighbor in adjacency[current]:
            if neighbor not in visited:
                queue.append((neighbor, path + [neighbor]))
    return []

# Generate paths for each box
paths = {box: find_path(initial_state[box], goal_state[box], adjacency) for box in initial_state}

# Initialize the plan with the initial state
plan = [initial_state.copy()]

# Function to move boxes according to their paths
def move_boxes(paths, plan):
    current_state = initial_state.copy()
    while current_state != goal_state:
        for box, path in paths.items():
            if current_state[box] != goal_state[box]:
                # Find the next step for the box
                current_index = path.index(current_state[box])
                next_position = path[current_index + 1]
                # Check if the next position is free
                if next_position not in current_state.values():
                    # Move the box
                    current_state[box] = next_position
                    plan.append(current_state.copy())
                    break

# Execute the moves
move_boxes(paths, plan)

# Print the plan
print(plan)