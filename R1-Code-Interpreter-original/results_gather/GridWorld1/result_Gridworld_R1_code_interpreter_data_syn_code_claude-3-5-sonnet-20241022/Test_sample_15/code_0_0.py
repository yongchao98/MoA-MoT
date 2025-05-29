from collections import deque
import json

def find_shortest_path(start, end, adjacency):
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
            
        for next_pos in adjacency[current]:
            if next_pos not in visited:
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_complete_path(start, goals, adjacency):
    current = start
    remaining_goals = set(goals)
    path = [start]
    
    while remaining_goals:
        # Find the closest remaining goal
        min_dist = float('inf')
        next_goal = None
        next_path = None
        
        for goal in remaining_goals:
            temp_path = find_shortest_path(current, goal, adjacency)
            if temp_path and len(temp_path) < min_dist:
                min_dist = len(temp_path)
                next_goal = goal
                next_path = temp_path
        
        if next_path is None:
            return None  # No valid path found
            
        # Add the path to this goal (excluding the start position)
        path.extend(next_path[1:])
        remaining_goals.remove(next_goal)
        current = next_goal
    
    return path

# Problem setup
start = "C3,1"
goals = ["C2,2", "C1,3", "C4,3", "C3,4"]
obstacles = ["C4,4", "C4,1", "C4,2", "C3,5", "C1,2"]

# Load adjacency data
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

# Remove obstacle connections
for obs in obstacles:
    if obs in adjacency:
        # Remove this position from its neighbors' adjacency lists
        for neighbor in adjacency[obs]:
            if neighbor in adjacency:
                adjacency[neighbor] = [pos for pos in adjacency[neighbor] if pos != obs]
        # Remove the obstacle position from the adjacency dict
        adjacency[obs] = []

# Find the complete path
result = find_complete_path(start, goals, adjacency)
print(f"<<<{json.dumps(result)}>>>")