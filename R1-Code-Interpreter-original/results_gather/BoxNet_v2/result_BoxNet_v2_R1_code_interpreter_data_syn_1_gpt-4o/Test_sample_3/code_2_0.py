import json
from collections import deque

# Initial and goal states
initial_state = {
    "box1": "C1,7",
    "box2": "C2,8",
    "box3": "C4,5",
    "box4": "C2,4",
    "box5": "C4,7",
    "box6": "C3,6"
}

goal_state = {
    "box1": "C3,6",
    "box2": "C4,1",
    "box3": "C2,7",
    "box4": "C3,3",
    "box5": "C2,3",
    "box6": "C2,6"
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
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"],
    "C2,6": ["C2,5", "C2,7", "C1,6", "C3,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7", "C3,7"],
    "C2,8": ["C2,7", "C1,8", "C3,8"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C3,6": ["C3,5", "C3,7", "C2,6", "C4,6"],
    "C3,7": ["C3,6", "C3,8", "C2,7", "C4,7"],
    "C3,8": ["C3,7", "C2,8", "C4,8"],
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"],
    "C4,6": ["C4,5", "C4,7", "C3,6"],
    "C4,7": ["C4,6", "C4,8", "C3,7"],
    "C4,8": ["C4,7", "C3,8"]
}

# Function to find the shortest path using BFS
def bfs_shortest_path(start, goal, adjacency):
    queue = deque([(start, [start])])
    visited = set()

    while queue:
        current, path = queue.popleft()
        if current == goal:
            return path
        if current not in visited:
            visited.add(current)
            for neighbor in adjacency[current]:
                if neighbor not in visited:
                    queue.append((neighbor, path + [neighbor]))
    return []

# Plan generation
def generate_plan(initial_state, goal_state, adjacency):
    plan = [initial_state.copy()]
    current_state = initial_state.copy()

    # Find paths for each box
    paths = {}
    for box, start in initial_state.items():
        goal = goal_state[box]
        if start != goal:
            paths[box] = bfs_shortest_path(start, goal, adjacency)

    # Move boxes along their paths
    max_steps = max(len(path) for path in paths.values())
    for step in range(1, max_steps):
        next_state = current_state.copy()
        for box, path in paths.items():
            if step < len(path):
                next_position = path[step]
                # Ensure no two boxes occupy the same cell
                if next_position not in next_state.values():
                    next_state[box] = next_position
        plan.append(next_state)
        current_state = next_state

    return plan

# Generate the plan
plan = generate_plan(initial_state, goal_state, adjacency)

# Output the plan
print(json.dumps(plan, indent=2))