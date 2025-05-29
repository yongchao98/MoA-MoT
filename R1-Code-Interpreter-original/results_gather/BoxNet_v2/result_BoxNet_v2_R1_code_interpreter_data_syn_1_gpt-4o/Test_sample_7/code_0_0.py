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

# BFS to find the plan
def find_plan(initial_state, goal_state, adjacency):
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(tuple(initial_state.items()))

    while queue:
        current_state, path = queue.popleft()

        if current_state == goal_state:
            return path + [current_state]

        for box, current_cell in current_state.items():
            if current_cell == goal_state[box]:
                continue

            for neighbor in adjacency[current_cell]:
                if neighbor not in current_state.values():
                    new_state = current_state.copy()
                    new_state[box] = neighbor
                    state_tuple = tuple(new_state.items())

                    if state_tuple not in visited:
                        visited.add(state_tuple)
                        queue.append((new_state, path + [current_state]))

    return []

plan = find_plan(initial_state, goal_state, adjacency)
print(plan)