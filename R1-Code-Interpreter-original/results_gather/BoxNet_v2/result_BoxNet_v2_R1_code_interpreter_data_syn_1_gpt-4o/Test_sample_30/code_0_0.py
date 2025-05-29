from collections import deque

# Initial and goal states
initial_state = {
    "box1": "C1,2",
    "box2": "C2,4",
    "box3": "C2,3",
    "box4": "C1,4"
}

goal_state = {
    "box1": "C3,1",
    "box2": "C2,3",
    "box3": "C3,2",
    "box4": "C1,5"
}

# Adjacency list
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
    "C3,1": ["C3,2", "C2,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5"],
    "C3,6": ["C3,5", "C2,6"]
}

# Function to check if a state is the goal state
def is_goal(state):
    return state == goal_state

# Function to generate successors of a given state
def generate_successors(state):
    successors = []
    for box, location in state.items():
        if location != goal_state[box]:  # Only move if not at goal
            for neighbor in adjacency[location]:
                # Check if the move is valid (no two boxes in the same cell)
                if neighbor not in state.values():
                    new_state = state.copy()
                    new_state[box] = neighbor
                    successors.append(new_state)
    return successors

# BFS to find the plan
def bfs(initial_state):
    queue = deque([(initial_state, [initial_state])])
    visited = set()
    visited.add(tuple(initial_state.items()))

    while queue:
        current_state, path = queue.popleft()

        if is_goal(current_state):
            return path

        for successor in generate_successors(current_state):
            successor_tuple = tuple(successor.items())
            if successor_tuple not in visited:
                visited.add(successor_tuple)
                queue.append((successor, path + [successor]))

    return None

# Generate the plan
plan = bfs(initial_state)
print(plan)