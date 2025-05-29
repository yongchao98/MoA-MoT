def generate_plan(initial_state, goal_state, adjacency):
    # Initialize the plan with the initial state
    plan = [initial_state.copy()]
    
    # Function to find a path for a box from start to goal
    def find_path(start, goal):
        from collections import deque
        
        # BFS to find the shortest path
        queue = deque([(start, [start])])
        visited = set()
        
        while queue:
            current, path = queue.popleft()
            if current == goal:
                return path
            
            for neighbor in adjacency[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, path + [neighbor]))
        
        return []
    
    # Move each box to its goal
    for box, start in initial_state.items():
        goal = goal_state[box]
        if start == goal:
            continue
        
        path = find_path(start, goal)
        
        # Move the box along the path
        for step in path[1:]:
            # Create a new state
            new_state = plan[-1].copy()
            new_state[box] = step
            plan.append(new_state)
    
    return plan

# Define the initial and goal states
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

# Define the adjacency list
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

# Generate the plan
plan = generate_plan(initial_state, goal_state, adjacency)

# Print the plan
print(plan)