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

# Helper function to calculate Manhattan distance
def manhattan_distance(cell1, cell2):
    r1, c1 = map(int, cell1[1:].split(','))
    r2, c2 = map(int, cell2[1:].split(','))
    return abs(r1 - r2) + abs(c1 - c2)

# A* search to find the path
def a_star_search(initial_state, goal_state, adjacency):
    queue = deque([(initial_state, [], 0)])
    visited = set()
    visited.add(tuple(initial_state.items()))

    while queue:
        current_state, path, cost = queue.popleft()

        if current_state == goal_state:
            return path + [current_state]

        for box, current_cell in current_state.items():
            if current_cell == goal_state[box]:
                continue

            for neighbor in adjacency[current_cell]:
                new_state = current_state.copy()
                new_state[box] = neighbor

                if len(set(new_state.values())) == len(new_state):  # Ensure no two boxes in the same cell
                    state_tuple = tuple(new_state.items())
                    if state_tuple not in visited:
                        visited.add(state_tuple)
                        new_cost = cost + 1
                        heuristic = sum(manhattan_distance(new_state[b], goal_state[b]) for b in new_state)
                        queue.append((new_state, path + [current_state], new_cost + heuristic))
                        queue = deque(sorted(queue, key=lambda x: x[2]))  # Sort by cost + heuristic

    return []

# Generate the plan
plan = a_star_search(initial_state, goal_state, adjacency)
print(plan)