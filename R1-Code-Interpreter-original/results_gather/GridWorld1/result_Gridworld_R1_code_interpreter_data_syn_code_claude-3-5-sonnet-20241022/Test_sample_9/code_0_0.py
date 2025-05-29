import json
from collections import deque
from itertools import permutations

def bfs_shortest_path(graph, start, end, obstacles):
    if start == end:
        return [start]
    
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        current = path[-1]
        
        for neighbor in graph[current]:
            if neighbor in obstacles:
                continue
            if neighbor == end:
                return path + [neighbor]
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(path + [neighbor])
    return None

def find_path_through_goals(graph, start, goals, obstacles):
    # Get all valid paths between each pair of points
    paths = {}
    all_points = [start] + goals
    
    for p1 in all_points:
        paths[p1] = {}
        for p2 in all_points:
            if p1 != p2:
                path = bfs_shortest_path(graph, p1, p2, obstacles)
                if path:
                    paths[p1][p2] = path

    # Try different goal orderings
    best_path = None
    min_length = float('inf')
    
    # Try all permutations of goals
    for goal_order in permutations(goals):
        current = start
        full_path = [start]
        valid_sequence = True
        
        for next_goal in goal_order:
            if next_goal not in paths[current]:
                valid_sequence = False
                break
            path_segment = paths[current][next_goal][1:]  # Exclude start point to avoid duplicates
            full_path.extend(path_segment)
            current = next_goal
            
        if valid_sequence and len(full_path) < min_length:
            min_length = len(full_path)
            best_path = full_path

    return best_path

# Main execution
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
    "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C3,5", "C5,5"],
    "C5,1": ["C5,2", "C4,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C4,5"]
}

start = "C1,3"
goals = ["C5,1", "C4,3", "C5,4", "C1,4", "C3,4", "C3,2"]
obstacles = ["C5,2", "C5,3", "C4,5", "C5,5", "C3,5", "C2,4", "C1,5", "C2,5"]

# Find the optimal path
path = find_path_through_goals(adjacency, start, goals, obstacles)

# Format the output
print(f"<<<{json.dumps(path)}>>>")