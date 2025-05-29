import json
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

# Plan generation using BFS
def generate_plan(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(tuple(initial_state.items()))

    while queue:
        current_state, path = queue.popleft()

        if current_state == goal_state:
            return path + [current_state]

        for box in initial_state.keys():
            if current_state[box] != goal_state[box]:
                for adj in adjacency[current_state[box]]:
                    if adj not in current_state.values() or adj == goal_state[box]:
                        new_state = current_state.copy()
                        new_state[box] = adj
                        state_tuple = tuple(new_state.items())
                        if state_tuple not in visited:
                            visited.add(state_tuple)
                            queue.append((new_state, path + [current_state]))

    return []

# Generate the plan
plan = generate_plan(initial_state, goal_state)

# Output the plan
print(json.dumps(plan, indent=2))