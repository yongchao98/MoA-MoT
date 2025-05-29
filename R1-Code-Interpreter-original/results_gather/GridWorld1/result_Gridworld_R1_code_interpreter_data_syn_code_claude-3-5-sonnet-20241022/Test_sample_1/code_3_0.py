import json
from collections import deque

def find_optimal_complete_path(start, goals, obstacles, adjacency):
    def bfs(start, target):
        queue = deque([(start, [start])])
        visited = {start}
        
        while queue:
            current, path = queue.popleft()
            if current == target:
                return path
            
            for next_pos in adjacency[current]:
                if next_pos not in visited and next_pos not in obstacles:
                    visited.add(next_pos)
                    queue.append((next_pos, path + [next_pos]))
        return None

    def optimize_path(path):
        # Remove cycles while preserving goal visits
        visited = set()
        goals_visited = set()
        optimized = []
        
        for pos in path:
            if pos not in visited or (pos in goals and pos not in goals_visited):
                optimized.append(pos)
                visited.add(pos)
                if pos in goals:
                    goals_visited.add(pos)
        return optimized

    # Find path through all goals
    current = start
    final_path = [start]
    remaining_goals = set(goals)
    visited_positions = {start}

    while remaining_goals:
        # Find nearest unvisited goal
        best_dist = float('inf')
        best_goal = None
        best_path = None
        
        for goal in remaining_goals:
            path = bfs(current, goal)
            if path and len(path) < best_dist:
                best_dist = len(path)
                best_goal = goal
                best_path = path
        
        if best_goal is None:
            return None
            
        # Add path to the nearest goal
        path_segment = best_path[1:]  # Skip the first position as it's already included
        final_path.extend(path_segment)
        
        # Update state
        current = best_goal
        remaining_goals.remove(best_goal)
        visited_positions.update(path_segment)

    # Optimize the final path
    return optimize_path(final_path)

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
path = find_optimal_complete_path(start, goals, obstacles, adjacency)

# Output the result
print(f"<<<{json.dumps(path)}>>>")