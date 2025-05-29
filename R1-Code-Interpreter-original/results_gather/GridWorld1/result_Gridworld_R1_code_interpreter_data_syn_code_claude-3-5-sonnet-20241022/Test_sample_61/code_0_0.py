import json
from collections import deque, defaultdict
from copy import deepcopy

def find_shortest_path(adjacency, start, end, obstacles):
    if start == end:
        return [start]
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                if next_pos == end:
                    return path + [next_pos]
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_complete_path(adjacency, start, goals, obstacles):
    remaining_goals = set(goals)
    current_pos = start
    final_path = [start]
    
    # If start position is a goal, mark it as visited
    if start in remaining_goals:
        remaining_goals.remove(start)
    
    while remaining_goals:
        # Find the nearest unvisited goal
        min_path = None
        min_goal = None
        min_length = float('inf')
        
        for goal in remaining_goals:
            path = find_shortest_path(adjacency, current_pos, goal, obstacles)
            if path and (min_path is None or len(path) < min_length):
                min_path = path
                min_goal = goal
                min_length = len(path)
        
        if min_path is None:
            return None  # No valid path found
        
        # Add the path to the final path (excluding the first position as it's already included)
        final_path.extend(min_path[1:])
        current_pos = min_goal
        remaining_goals.remove(min_goal)
    
    return final_path

# Main execution
start = "C5,6"
goals = ['C5,3', 'C5,2', 'C4,1', 'C2,1', 'C4,4', 'C6,5', 'C2,3']
obstacles = ['C3,5', 'C6,2', 'C3,4', 'C1,2', 'C1,4', 'C3,3', 'C4,5', 'C5,4', 'C3,2', 'C1,3']

# Parse the adjacency dictionary
adjacency = {
    # ... (using the provided adjacency dictionary)
}

# Find the complete path
path = find_complete_path(adjacency, start, goals, obstacles)

# Format and print the result
if path:
    result = json.dumps(path)
    print(f"<<<{result}>>>")
else:
    print("No valid path found")