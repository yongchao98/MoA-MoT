import heapq

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

# Convert cell to coordinates
def cell_to_coords(cell):
    row, col = map(int, cell[1:].split(','))
    return row, col

# Heuristic function: Manhattan distance
def heuristic(state, goal_state):
    return sum(abs(cell_to_coords(state[box])[0] - cell_to_coords(goal_state[box])[0]) +
               abs(cell_to_coords(state[box])[1] - cell_to_coords(goal_state[box])[1])
               for box in state)

# A* algorithm to generate the plan
def a_star_plan(initial_state, goal_state, adjacency):
    open_set = []
    heapq.heappush(open_set, (0, tuple(sorted(initial_state.items())), []))
    visited = set()

    while open_set:
        _, current_state_tuple, path = heapq.heappop(open_set)
        current_state = dict(current_state_tuple)

        if current_state_tuple in visited:
            continue
        visited.add(current_state_tuple)

        # Check if the current state is the goal state
        if current_state == goal_state:
            return path + [current_state]

        # Generate new states by moving boxes
        for box, current_cell in current_state.items():
            if current_cell == goal_state[box]:
                continue  # Skip if the box is already at its goal

            for adjacent in adjacency[current_cell]:
                # Check if the adjacent cell is free
                if adjacent not in current_state.values():
                    new_state = current_state.copy()
                    new_state[box] = adjacent
                    new_cost = len(path) + 1 + heuristic(new_state, goal_state)
                    heapq.heappush(open_set, (new_cost, tuple(sorted(new_state.items())), path + [current_state]))

    return []

# Generate the plan
plan = a_star_plan(initial_state, goal_state, adjacency)
print(plan)