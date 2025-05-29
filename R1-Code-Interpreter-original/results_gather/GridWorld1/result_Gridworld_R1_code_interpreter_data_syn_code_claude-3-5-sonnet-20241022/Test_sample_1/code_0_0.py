import json
from collections import deque

def find_path_visiting_all_goals(start, goals, obstacles, adjacency):
    def bfs(start, target, visited_positions):
        queue = deque([(start, [start])])
        visited = set([start])
        
        while queue:
            current, path = queue.popleft()
            if current == target:
                return path
                
            for next_pos in adjacency[current]:
                if next_pos not in visited and next_pos not in obstacles:
                    visited.add(next_pos)
                    new_path = path + [next_pos]
                    queue.append((next_pos, new_path))
        return None

    remaining_goals = set(goals)
    current_pos = start
    final_path = [start]
    visited_goals = set()

    while remaining_goals:
        # Find the nearest unvisited goal
        min_path = None
        nearest_goal = None
        
        for goal in remaining_goals:
            path = bfs(current_pos, goal, visited_goals)
            if path and (min_path is None or len(path) < len(min_path)):
                min_path = path
                nearest_goal = goal
        
        if min_path is None:
            return None  # No valid path found
            
        # Add the path to the nearest goal (excluding the start position)
        final_path.extend(min_path[1:])
        current_pos = nearest_goal
        remaining_goals.remove(nearest_goal)
        visited_goals.add(nearest_goal)

    return final_path

# Define the problem parameters
start = "C2,5"
goals = ["C3,1", "C4,4", "C4,5", "C1,4"]
obstacles = ["C3,4", "C1,5", "C1,2", "C2,2", "C4,1"]
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

# Find the path
path = find_path_visiting_all_goals(start, goals, obstacles, adjacency)

# Output the result
print(f"<<<{json.dumps(path)}>>>")