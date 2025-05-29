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

# Function to generate the plan
def generate_plan(initial_state, goal_state, adjacency):
    plan = [initial_state.copy()]
    current_state = initial_state.copy()
    
    # Function to move a single box to its goal
    def move_box(box, current_location, goal_location):
        queue = deque([(current_location, [current_location])])
        visited = set()
        
        while queue:
            position, path = queue.popleft()
            if position == goal_location:
                return path
            
            for neighbor in adjacency[position]:
                if neighbor not in visited and neighbor not in current_state.values():
                    visited.add(neighbor)
                    queue.append((neighbor, path + [neighbor]))
        
        return [current_location]  # No move possible
    
    # Move each box to its goal
    for box in initial_state.keys():
        path = move_box(box, current_state[box], goal_state[box])
        for step in path[1:]:
            current_state[box] = step
            plan.append(current_state.copy())
    
    return plan

# Generate the plan
plan = generate_plan(initial_state, goal_state, adjacency)
print(plan)