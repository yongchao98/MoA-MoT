import json
from collections import deque

def get_neighbors(pos, adjacency, obstacles, visited):
    return [n for n in adjacency[pos] if n not in obstacles and n not in visited]

def find_path_to_closest_goal(start, goals, adjacency, obstacles, visited):
    queue = deque([(start, [start])])
    seen = {start}
    
    while queue:
        current, path = queue.popleft()
        
        if current in goals:
            return path
        
        for next_pos in adjacency[current]:
            if next_pos not in obstacles and next_pos not in seen:
                seen.add(next_pos)
                new_path = path + [next_pos]
                queue.append((next_pos, new_path))
    
    return None

def optimize_path(path, goals, adjacency, obstacles):
    result = [path[0]]
    i = 0
    while i < len(path) - 1:
        current = path[i]
        # Look ahead for the next goal
        next_goal_idx = i + 1
        while next_goal_idx < len(path) and path[next_goal_idx] not in goals:
            next_goal_idx += 1
            
        if next_goal_idx < len(path):
            # Try to find a direct path to the next goal
            direct_path = find_path_to_closest_goal(
                current,
                {path[next_goal_idx]},
                adjacency,
                obstacles,
                set(result[:-1])  # Avoid previously visited positions except the current
            )
            
            if direct_path:
                result.extend(direct_path[1:])
                i = next_goal_idx
            else:
                result.append(path[i + 1])
                i += 1
        else:
            result.append(path[i + 1])
            i += 1
            
    return result

def find_optimal_route(start, goals, adjacency, obstacles):
    unvisited_goals = set(goals)
    current = start
    path = [start]
    visited = {start}
    
    while unvisited_goals:
        # Find the closest unvisited goal
        best_path = None
        best_goal = None
        min_length = float('inf')
        
        for goal in unvisited_goals:
            temp_path = find_path_to_closest_goal(
                current,
                {goal},
                adjacency,
                obstacles,
                visited - {current}  # Allow revisiting current position
            )
            
            if temp_path and len(temp_path) < min_length:
                min_length = len(temp_path)
                best_path = temp_path
                best_goal = goal
        
        if not best_path:
            return None
        
        # Add the path to the goal
        path.extend(best_path[1:])
        current = best_goal
        visited.update(best_path)
        unvisited_goals.remove(best_goal)
    
    # Optimize the final path
    return optimize_path(path, set(goals), adjacency, obstacles)

# Problem setup
start = "C4,4"
goals = ['C4,2', 'C5,1', 'C5,4', 'C3,2', 'C3,3', 'C5,3']
obstacles = ['C1,2', 'C2,5', 'C5,2', 'C1,5', 'C3,5', 'C1,3', 'C1,4', 'C4,3']

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

solution = find_optimal_route(start, goals, adjacency, obstacles)
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No valid path found")