import json
from collections import deque, defaultdict
import heapq

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

def find_all_paths_between_goals(start, goals, adjacency, obstacles):
    all_points = goals + [start]
    paths = {}
    
    for i, point1 in enumerate(all_points):
        for point2 in all_points[i+1:]:
            path = find_shortest_path(point1, point2, adjacency, obstacles)
            if path:
                paths[(point1, point2)] = path
                paths[(point2, point1)] = path[::-1]
    return paths

def find_solution_path(start, goals, adjacency, obstacles):
    # Get all shortest paths between points
    paths = find_all_paths_between_goals(start, goals, adjacency, obstacles)
    
    # Try to find a good order to visit goals
    unvisited = set(goals)
    current = start
    final_path = [start]
    
    while unvisited:
        # Find the closest unvisited goal
        best_dist = float('inf')
        best_goal = None
        best_path = None
        
        for goal in unvisited:
            if (current, goal) in paths:
                path = paths[(current, goal)]
                if len(path) < best_dist:
                    best_dist = len(path)
                    best_goal = goal
                    best_path = path
        
        if best_goal is None:
            return None  # No valid solution found
            
        # Add the path to this goal (excluding the start position)
        final_path.extend(best_path[1:])
        current = best_goal
        unvisited.remove(best_goal)
    
    return final_path

# Main execution
start = "C2,5"
goals = ['C1,4', 'C3,2', 'C1,5', 'C4,4', 'C1,2', 'C2,2']
obstacles = ['C4,3', 'C3,5', 'C1,1', 'C2,3', 'C5,3', 'C2,1', 'C4,2', 'C5,1']

# Parse adjacency from the given dictionary
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

solution = find_solution_path(start, goals, adjacency, obstacles)
print(f"<<<{json.dumps(solution)}>>>")