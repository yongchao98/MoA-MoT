from collections import deque
import json

def find_shortest_path(start, end, adjacency, obstacles):
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

def find_complete_path():
    # Initialize data
    start = "C1,2"
    goals = ['C4,1', 'C3,3', 'C5,1', 'C2,3', 'C1,4', 'C5,2', 'C3,2']
    obstacles = ['C1,6', 'C3,6', 'C4,3', 'C5,3', 'C5,6', 'C2,6', 'C4,4', 'C3,5', 'C4,2', 'C4,5']
    
    adjacency = {
        # ... (using the provided adjacency dictionary)
    }
    
    # Initialize variables
    current_pos = start
    remaining_goals = set(goals)
    final_path = [start]
    
    # Keep finding paths to nearest goals until all goals are visited
    while remaining_goals:
        # Find the nearest unvisited goal
        shortest_path = None
        nearest_goal = None
        shortest_length = float('inf')
        
        for goal in remaining_goals:
            path = find_shortest_path(current_pos, goal, adjacency, obstacles)
            if path and len(path) < shortest_length:
                shortest_length = len(path)
                shortest_path = path
                nearest_goal = goal
        
        if not shortest_path:
            return None  # No valid path found
        
        # Add the path to the final path (excluding the start position)
        final_path.extend(shortest_path[1:])
        current_pos = nearest_goal
        remaining_goals.remove(nearest_goal)
    
    return final_path

# Find and print the path
path = find_complete_path()
if path:
    print(f"<<<{json.dumps(path)}>>>")
else:
    print("No valid path found")