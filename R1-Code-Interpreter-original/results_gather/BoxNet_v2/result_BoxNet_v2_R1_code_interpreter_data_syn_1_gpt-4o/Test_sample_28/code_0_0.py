from collections import deque

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

# Function to find path for a single box
def find_path(start, goal, occupied):
    queue = deque([(start, [start])])
    visited = set()
    while queue:
        current, path = queue.popleft()
        if current == goal:
            return path
        for neighbor in adjacency[current]:
            if neighbor not in visited and neighbor not in occupied:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    return []

# Generate the plan
def generate_plan(initial, goal):
    plan = [initial.copy()]
    current_state = initial.copy()
    while current_state != goal:
        for box, current_pos in current_state.items():
            if current_pos != goal[box]:
                path = find_path(current_pos, goal[box], set(current_state.values()))
                if len(path) > 1:
                    current_state[box] = path[1]
                    plan.append(current_state.copy())
                    break
    return plan

# Generate and print the plan
plan = generate_plan(initial_state, goal_state)
print(plan)