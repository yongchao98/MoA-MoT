import json
from collections import deque, defaultdict

def find_path(start, goals, obstacles, adjacency):
    def bfs(start, target):
        queue = deque([(start, [start])])
        visited = {start}
        
        while queue:
            current, path = queue.popleft()
            if current == target:
                return path
                
            for next_cell in adjacency[current]:
                if next_cell not in visited and next_cell not in obstacles:
                    visited.add(next_cell)
                    queue.append((next_cell, path + [next_cell]))
        return None

    # Initialize variables
    current_pos = start
    final_path = [start]
    unvisited_goals = set(goals)
    
    # Keep finding paths to nearest unvisited goals until all are visited
    while unvisited_goals:
        # Find distances to all remaining goals
        min_dist = float('inf')
        next_goal = None
        next_path = None
        
        for goal in unvisited_goals:
            path = bfs(current_pos, goal)
            if path and len(path) < min_dist:
                min_dist = len(path)
                next_goal = goal
                next_path = path
        
        if next_path is None:
            return None  # No valid path found
            
        # Add the path to the next goal (excluding the start position)
        final_path.extend(next_path[1:])
        current_pos = next_goal
        unvisited_goals.remove(next_goal)
    
    print(json.dumps(final_path))

# Problem setup
start = "C4,2"
goals = ["C2,2", "C3,2", "C1,2", "C3,1"]
obstacles = ["C1,1", "C3,4", "C2,4", "C4,4", "C2,3"]
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C2,4"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C1,4", "C3,4"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C2,4", "C4,4"],
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C3,4"]
}

find_path(start, goals, obstacles, adjacency)