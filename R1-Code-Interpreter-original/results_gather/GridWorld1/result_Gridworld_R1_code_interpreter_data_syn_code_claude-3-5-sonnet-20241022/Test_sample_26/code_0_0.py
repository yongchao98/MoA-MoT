from collections import deque
import json

def find_shortest_path(adjacency, start, end, obstacles):
    if start == end:
        return [start]
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                if next_pos == end:
                    return path + [next_pos]
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_complete_path(adjacency, start, goals, obstacles):
    current_pos = start
    remaining_goals = set(goals)
    path = [start]
    
    while remaining_goals:
        # Find the nearest unvisited goal
        min_distance = float('inf')
        next_goal = None
        next_path = None
        
        for goal in remaining_goals:
            current_path = find_shortest_path(adjacency, current_pos, goal, obstacles)
            if current_path and len(current_path) < min_distance:
                min_distance = len(current_path)
                next_goal = goal
                next_path = current_path
        
        if next_goal is None:
            return None  # No valid path found
            
        # Add the path to the next goal (excluding the start position)
        path.extend(next_path[1:])
        current_pos = next_goal
        remaining_goals.remove(next_goal)
    
    return path

# Problem setup
start = "C2,6"
goals = ['C6,3', 'C4,2', 'C2,5', 'C4,3', 'C4,5', 'C1,1', 'C3,3']
obstacles = ['C5,6', 'C4,6', 'C1,4', 'C1,3', 'C6,2', 'C3,4', 'C6,1', 'C3,2', 'C4,1', 'C6,6']

# Load adjacency data
adjacency = {
    # ... (your provided adjacency dictionary)
}

# Find the solution path
solution_path = find_complete_path(adjacency, start, goals, obstacles)

# Format the output
if solution_path:
    print(f"<<<{json.dumps(solution_path)}>>>")
else:
    print("No valid path found")