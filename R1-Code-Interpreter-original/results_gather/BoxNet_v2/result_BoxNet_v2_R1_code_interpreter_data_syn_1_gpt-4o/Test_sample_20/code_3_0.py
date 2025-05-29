from collections import deque

# Initial and goal states
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

# Adjacency list
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

# Function to find the shortest path using BFS
def bfs(start, goal):
    queue = deque([(start, [start])])
    visited = set()
    while queue:
        current, path = queue.popleft()
        if current == goal:
            return path
        if current in visited:
            continue
        visited.add(current)
        for neighbor in adjacency[current]:
            if neighbor not in visited:
                queue.append((neighbor, path + [neighbor]))
    return []

# Generate paths for each box
paths = {box: bfs(initial_state[box], goal_state[box]) for box in initial_state}

# Generate the plan
plan = [initial_state.copy()]
current_state = initial_state.copy()

# Move boxes along their paths
while current_state != goal_state:
    for box, path in paths.items():
        if current_state[box] != goal_state[box]:
            # Move the box one step along its path
            next_position = path[path.index(current_state[box]) + 1]
            current_state[box] = next_position
            plan.append(current_state.copy())
            break

# Ensure the final state matches the goal state
assert current_state == goal_state, "Final state does not match goal state"

# Print the plan in a controlled manner
for state in plan:
    print(state)