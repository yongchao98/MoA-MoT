import json
from collections import deque

# Define the problem parameters
initial_position = "C1,2"
goals = {'C1,4', 'C4,4', 'C2,1', 'C3,1', 'C3,4', 'C3,5'}
obstacles = {'C4,2', 'C2,2', 'C3,2', 'C5,2', 'C4,3', 'C3,3', 'C5,3', 'C2,3'}
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
    "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C3,5", "C5,5"],
    "C5,1": ["C5,2", "C4,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C4,5"]
}

# Initialize the queue with the initial position
queue = deque([(initial_position, [initial_position], set())])

# Perform a breadth-first search
while queue:
    current_position, path, visited_goals = queue.popleft()
    
    # Check if all goals have been visited
    if visited_goals == goals:
        print(json.dumps(path))
        break
    
    # Explore adjacent cells
    for neighbor in adjacency.get(current_position, []):
        if neighbor not in obstacles:
            new_visited_goals = visited_goals.copy()
            if neighbor in goals:
                new_visited_goals.add(neighbor)
            queue.append((neighbor, path + [neighbor], new_visited_goals))