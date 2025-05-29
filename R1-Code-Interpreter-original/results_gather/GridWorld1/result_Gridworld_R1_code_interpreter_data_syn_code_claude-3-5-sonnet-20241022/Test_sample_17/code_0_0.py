import json
from collections import deque
from itertools import permutations

def find_path(start, goals, obstacles, adjacency):
    def bfs(start, target):
        if start == target:
            return [start]
        visited = {start}
        queue = deque([(start, [start])])
        
        while queue:
            current, path = queue.popleft()
            for next_pos in adjacency[current]:
                if next_pos not in visited and next_pos not in obstacles:
                    if next_pos == target:
                        return path + [next_pos]
                    visited.add(next_pos)
                    queue.append((next_pos, path + [next_pos]))
        return None

    # Find paths between all pairs of relevant points
    points = [start] + goals
    paths = {}
    for p1 in points:
        for p2 in points:
            if p1 != p2:
                path = bfs(p1, p2)
                if path:
                    paths[(p1, p2)] = path

    # Try different goal orderings
    best_path = None
    min_length = float('inf')
    
    # Use a greedy approach instead of trying all permutations
    current = start
    remaining_goals = set(goals)
    final_path = [start]
    
    while remaining_goals:
        # Find the closest unvisited goal
        best_next = None
        best_next_path = None
        best_dist = float('inf')
        
        for goal in remaining_goals:
            if (current, goal) in paths:
                path = paths[(current, goal)]
                if len(path) < best_dist:
                    best_dist = len(path)
                    best_next = goal
                    best_next_path = path
        
        if best_next is None:
            return None  # No valid path found
            
        # Add the path to the next goal (excluding the start point to avoid duplication)
        final_path.extend(best_next_path[1:])
        current = best_next
        remaining_goals.remove(best_next)

    return final_path

# Main execution
start = "C3,4"
goals = ['C3,1', 'C5,2', 'C4,4', 'C1,2', 'C5,1', 'C4,3']
obstacles = ['C2,5', 'C3,5', 'C3,2', 'C2,2', 'C4,2', 'C2,4', 'C5,5', 'C4,5']

# Parse adjacency from the provided dictionary
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

path = find_path(start, goals, obstacles, adjacency)
print(f"<<<{json.dumps(path)}>>>")