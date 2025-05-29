from collections import deque

# Define the initial and goal states
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

# Define the adjacency list
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

# Function to generate the plan using BFS
def generate_plan_bfs(initial_state, goal_state, adjacency):
    queue = deque([(initial_state, [])])
    visited = set()

    while queue:
        current_state, path = queue.popleft()
        state_tuple = tuple(sorted(current_state.items()))

        if state_tuple in visited:
            continue

        visited.add(state_tuple)

        if current_state == goal_state:
            return path + [current_state]

        for box, current_location in current_state.items():
            if current_location != goal_state[box]:
                for neighbor in adjacency[current_location]:
                    new_state = current_state.copy()
                    new_state[box] = neighbor
                    queue.append((new_state, path + [current_state]))

    return []

# Generate the plan
plan = generate_plan_bfs(initial_state, goal_state, adjacency)
print(plan)