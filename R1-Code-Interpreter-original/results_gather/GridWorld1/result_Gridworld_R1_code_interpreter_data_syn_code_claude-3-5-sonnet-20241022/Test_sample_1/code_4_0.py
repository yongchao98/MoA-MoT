import json
from collections import deque
from itertools import permutations

def find_shortest_complete_path(start, goals, obstacles, adjacency):
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

    def get_path_through_goals(goal_order):
        current = start
        complete_path = [start]
        
        for goal in goal_order:
            path = bfs(current, goal)
            if not path:
                return None, float('inf')
            complete_path.extend(path[1:])  # Skip the first position as it's already included
            current = goal
            
        return complete_path, len(complete_path)

    # Try different goal orderings to find the shortest valid path
    best_path = None
    best_length = float('inf')
    
    # First try nearest neighbor approach
    current = start
    remaining_goals = goals.copy()
    ordered_goals = []
    
    while remaining_goals:
        nearest = None
        shortest_dist = float('inf')
        shortest_path = None
        
        for goal in remaining_goals:
            path = bfs(current, goal)
            if path and len(path) < shortest_dist:
                shortest_dist = len(path)
                nearest = goal
                shortest_path = path
        
        if nearest is None:
            break
            
        ordered_goals.append(nearest)
        remaining_goals.remove(nearest)
        current = nearest

    # Try the nearest neighbor solution
    path, length = get_path_through_goals(ordered_goals)
    if path and length < best_length:
        best_path = path
        best_length = length

    # Try some alternative permutations of goals
    for goal_order in permutations(goals, len(goals)):
        path, length = get_path_through_goals(goal_order)
        if path and length < best_length:
            best_path = path
            best_length = length

    return best_path

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
path = find_shortest_complete_path(start, goals, obstacles, adjacency)

# Output the result
print(f"<<<{json.dumps(path)}>>>")