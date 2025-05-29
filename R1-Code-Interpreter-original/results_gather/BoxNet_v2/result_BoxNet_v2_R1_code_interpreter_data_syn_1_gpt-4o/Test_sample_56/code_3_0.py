from collections import deque
import heapq

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

def bfs_shortest_path(start, goal, adjacency):
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

def heuristic(state, goal_state):
    return sum(1 for box, location in state.items() if location != goal_state[box])

def generate_plan(initial_state, goal_state, adjacency):
    plan = [initial_state.copy()]
    current_state = initial_state.copy()
    visited_states = set()

    while current_state != goal_state:
        for box, current_location in current_state.items():
            if current_location != goal_state[box]:
                path = bfs_shortest_path(current_location, goal_state[box], adjacency)
                if len(path) > 1:
                    next_location = path[1]
                    if next_location not in current_state.values():
                        current_state[box] = next_location
                        state_tuple = tuple(sorted(current_state.items()))
                        if state_tuple not in visited_states:
                            visited_states.add(state_tuple)
                            plan.append(current_state.copy())
                            break
    return plan

plan = generate_plan(initial_state, goal_state, adjacency)
print(plan)