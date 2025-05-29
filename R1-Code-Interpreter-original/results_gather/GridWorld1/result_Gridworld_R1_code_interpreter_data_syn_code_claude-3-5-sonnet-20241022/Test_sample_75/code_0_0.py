from collections import deque
import json

def find_shortest_path(start, end, adjacency, obstacles):
    if start in obstacles or end in obstacles:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
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
    current_pos = start
    remaining_goals = set(goals)
    path = [start]
    visited_goals = set()
    
    while remaining_goals:
        # Find the nearest unvisited goal
        best_path = None
        best_goal = None
        for goal in remaining_goals:
            current_path = find_shortest_path(current_pos, goal, adjacency, obstacles)
            if current_path is not None:
                if best_path is None or len(current_path) < len(best_path):
                    best_path = current_path
                    best_goal = goal
        
        if best_path is None:
            return None  # No valid path found
            
        # Add the path to the goal (excluding the start position if not the first goal)
        path.extend(best_path[1:])
        current_pos = best_goal
        remaining_goals.remove(best_goal)
        visited_goals.add(best_goal)
    
    return path

# Define the problem parameters
start = "C5,4"
goals = ['C6,5', 'C2,6', 'C6,6', 'C1,6', 'C1,4', 'C6,1', 'C6,2']
obstacles = ['C3,3', 'C4,6', 'C2,5', 'C4,5', 'C6,4', 'C3,5', 'C2,2', 'C3,4', 'C5,2', 'C2,4']

# Define adjacency (assuming it's provided in the full problem)
adjacency = {
    # ... (the full adjacency dictionary as provided in the problem)
}

# Find the complete path
complete_path = find_complete_path(start, goals, adjacency, obstacles)

# Format and print the result
if complete_path:
    result = json.dumps(complete_path)
    print(f"<<<{result}>>>")
else:
    print("No valid path found")