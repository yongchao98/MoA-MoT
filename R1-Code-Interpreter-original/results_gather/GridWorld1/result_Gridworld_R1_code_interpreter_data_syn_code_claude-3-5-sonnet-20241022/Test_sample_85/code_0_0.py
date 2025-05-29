import json
from collections import deque
from itertools import permutations

def find_shortest_path(adjacency, start, end, obstacles):
    if start == end:
        return [start]
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        for next_pos in adjacency[current]:
            if next_pos in obstacles:
                continue
            if next_pos == end:
                return path + [next_pos]
            if next_pos not in visited:
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_complete_path(adjacency, start, goals, obstacles):
    # Find all shortest paths between each pair of points
    points = [start] + goals
    paths = {}
    for p1 in points:
        for p2 in points:
            if p1 != p2:
                path = find_shortest_path(adjacency, p1, p2, obstacles)
                if path:
                    paths[(p1, p2)] = path

    # Try all possible goal orderings
    best_path = None
    min_length = float('inf')
    
    for goal_order in permutations(goals):
        # Check if we can reach all goals in this order
        current = start
        total_path = [start]
        valid_sequence = True
        
        for goal in goal_order:
            if (current, goal) not in paths:
                valid_sequence = False
                break
            path_segment = paths[(current, goal)][1:]  # Exclude start position
            total_path.extend(path_segment)
            current = goal
            
        if valid_sequence and len(total_path) < min_length:
            min_length = len(total_path)
            best_path = total_path

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
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C3,5"]
}

start = "C4,2"
goals = ["C2,5", "C1,3", "C2,1", "C3,3"]
obstacles = ["C3,2", "C1,1", "C3,5", "C3,1", "C1,2"]

path = find_complete_path(adjacency, start, goals, obstacles)
print(json.dumps(path))