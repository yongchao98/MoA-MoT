from collections import deque

# Define the problem parameters
initial_position = "C2,4"
goals = {'C2,3', 'C1,2', 'C2,5', 'C2,1'}
obstacles = {'C2,2', 'C4,3', 'C3,4', 'C1,4', 'C3,1'}
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C1,5", "C3,5"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"],
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C3,5"]
}

# Initialize the queue with the initial position
queue = deque([(initial_position, [initial_position], set())])

# Perform BFS to find a valid path
while queue:
    current_position, path, visited_goals = queue.popleft()
    
    # Check if all goals are visited
    if visited_goals == goals:
        print(path)
        break
    
    # Explore adjacent positions
    for neighbor in adjacency.get(current_position, []):
        if neighbor not in obstacles and neighbor not in path:
            new_visited_goals = visited_goals.copy()
            if neighbor in goals:
                new_visited_goals.add(neighbor)
            queue.append((neighbor, path + [neighbor], new_visited_goals))