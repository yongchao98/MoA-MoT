from collections import deque

# Define the initial state and goal state
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

# Function to check if a state is the goal state
def is_goal(state):
    return state == goal_state

# Function to generate the next possible states
def generate_next_states(current_state):
    next_states = []
    for box, current_pos in current_state.items():
        if current_pos != goal_state[box]:  # Only move if not at goal
            for neighbor in adjacency[current_pos]:
                # Check if the neighbor cell is free
                if neighbor not in current_state.values():
                    new_state = current_state.copy()
                    new_state[box] = neighbor
                    next_states.append(new_state)
    return next_states

# BFS to find the sequence of states
def bfs(initial_state):
    queue = deque([(initial_state, [initial_state])])
    visited = set()
    visited.add(tuple(initial_state.items()))

    while queue:
        current_state, path = queue.popleft()

        if is_goal(current_state):
            return path

        for next_state in generate_next_states(current_state):
            state_tuple = tuple(next_state.items())
            if state_tuple not in visited:
                visited.add(state_tuple)
                queue.append((next_state, path + [next_state]))

# Find the plan
plan = bfs(initial_state)

# Print the complete plan in the required format
print("<<<" + str(plan) + ">>>")