import json
from collections import deque

def find_optimal_path_no_backtrack(start, goals, obstacles, adjacency):
    def bfs_with_visited(start, target, blocked_positions=None):
        if blocked_positions is None:
            blocked_positions = set()
        queue = deque([(start, [start])])
        visited = {start}
        
        while queue:
            current, path = queue.popleft()
            if current == target:
                return path
            
            for next_pos in adjacency[current]:
                if (next_pos not in visited and 
                    next_pos not in obstacles and 
                    next_pos not in blocked_positions):
                    visited.add(next_pos)
                    queue.append((next_pos, path + [next_pos]))
        return None

    def find_path_to_all_goals():
        paths = []
        current = start
        remaining_goals = set(goals)
        visited = set([start])
        current_path = [start]
        
        while remaining_goals:
            best_path = None
            best_goal = None
            best_length = float('inf')
            
            # Try to find the nearest unvisited goal
            for goal in remaining_goals:
                path = bfs_with_visited(current, goal, visited)
                if path and len(path) < best_length:
                    best_length = len(path)
                    best_path = path
                    best_goal = goal
            
            if best_goal is None:
                return None
            
            # Add the path to the nearest goal (excluding start position)
            path_segment = best_path[1:]
            current_path.extend(path_segment)
            visited.update(path_segment)
            current = best_goal
            remaining_goals.remove(best_goal)
        
        return current_path

    # Get the complete path
    path = find_path_to_all_goals()
    
    if path:
        # Remove any remaining duplicates while preserving goal visits
        final_path = []
        visited = set()
        goals_visited = set()
        
        for pos in path:
            if pos not in visited or (pos in goals and pos not in goals_visited):
                final_path.append(pos)
                visited.add(pos)
                if pos in goals:
                    goals_visited.add(pos)
        
        return final_path
    return None

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
path = find_optimal_path_no_backtrack(start, goals, obstacles, adjacency)

# Output the result
print(f"<<<{json.dumps(path)}>>>")