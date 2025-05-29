from collections import deque
import json

def find_path(start, end, adjacency, obstacles):
    queue = deque([(start, [start])])
    visited = set([start])
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
            
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_complete_path(start, goals, adjacency, obstacles):
    current = start
    remaining_goals = set(goals)
    final_path = [start]
    
    while remaining_goals:
        # Find the nearest unvisited goal
        shortest_path = None
        nearest_goal = None
        shortest_length = float('inf')
        
        for goal in remaining_goals:
            path = find_path(current, goal, adjacency, obstacles)
            if path and len(path) < shortest_length:
                shortest_length = len(path)
                shortest_path = path
                nearest_goal = goal
        
        if not shortest_path:
            return None  # No valid path found
            
        # Add the path to the final path (excluding the start position)
        final_path.extend(shortest_path[1:])
        current = nearest_goal
        remaining_goals.remove(nearest_goal)
    
    return final_path

# Problem setup
start = "C3,3"
goals = ["C4,1", "C4,2", "C2,2", "C1,2"]
obstacles = ["C1,5", "C3,4", "C1,1", "C4,3", "C4,5"]

# Load adjacency data
adjacency = {
    # ... (your provided adjacency dict)
}

# Find the solution
solution = find_complete_path(start, goals, adjacency, obstacles)

# Output the solution
print(f"<<<{json.dumps(solution)}>>>")