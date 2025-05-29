import json
from collections import deque

def find_complete_path(start, goals, obstacles, adjacency):
    def bfs(start, target, visited_set=None):
        if visited_set is None:
            visited_set = set()
        queue = deque([(start, [start])])
        visited = {start}
        
        while queue:
            current, path = queue.popleft()
            if current == target:
                return path
            
            for next_pos in adjacency[current]:
                if (next_pos not in visited and 
                    next_pos not in obstacles and 
                    next_pos not in visited_set):
                    visited.add(next_pos)
                    queue.append((next_pos, path + [next_pos]))
        return None

    def find_best_goal_sequence():
        unvisited = goals.copy()
        current = start
        sequence = []
        while unvisited:
            best_dist = float('inf')
            best_goal = None
            best_path = None
            
            for goal in unvisited:
                path = bfs(current, goal)
                if path and len(path) < best_dist:
                    best_dist = len(path)
                    best_goal = goal
                    best_path = path
            
            if best_goal is None:
                return None
            
            sequence.append((best_goal, best_path))
            current = best_goal
            unvisited.remove(best_goal)
        return sequence

    # Find the complete path through all goals
    goal_sequence = find_best_goal_sequence()
    if not goal_sequence:
        return None

    # Construct the final path
    final_path = [start]
    for _, path in goal_sequence:
        final_path.extend(path[1:])  # Skip the first position as it's already included

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
path = find_complete_path(start, goals, obstacles, adjacency)

# Output the result
print(f"<<<{json.dumps(path)}>>>")