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

def optimize_path(path, goals):
    # Remove redundant moves while keeping all goals visited
    optimized = [path[0]]
    visited_goals = set()
    
    for i in range(1, len(path)):
        current = path[i]
        if current in goals:
            visited_goals.add(current)
        if i == len(path) - 1 or current in goals or path[i + 1] in goals:
            optimized.append(current)
    return optimized

def find_complete_path(adjacency, start, goals, obstacles):
    current_pos = start
    remaining_goals = set(goals)
    path = [start]
    visited_goals = set()
    
    while remaining_goals:
        min_cost = float('inf')
        next_goal = None
        next_path = None
        
        for goal in remaining_goals:
            current_path = find_shortest_path(adjacency, current_pos, goal, obstacles)
            if current_path:
                # Consider both distance and proximity to other goals
                cost = len(current_path)
                if cost < min_cost:
                    min_cost = cost
                    next_goal = goal
                    next_path = current_path
        
        if next_goal is None:
            return None
            
        path.extend(next_path[1:])
        current_pos = next_goal
        remaining_goals.remove(next_goal)
        visited_goals.add(next_goal)
    
    return optimize_path(path, goals)

# Problem setup
start = "C2,6"
goals = ['C6,3', 'C4,2', 'C2,5', 'C4,3', 'C4,5', 'C1,1', 'C3,3']
obstacles = ['C5,6', 'C4,6', 'C1,4', 'C1,3', 'C6,2', 'C3,4', 'C6,1', 'C3,2', 'C4,1', 'C6,6']

# Adjacency dictionary
adjacency = {
    # ... (same adjacency dictionary as before)
}

# Find the solution path
solution_path = find_complete_path(adjacency, start, goals, obstacles)

# Format the output
if solution_path:
    print(f"<<<{json.dumps(solution_path)}>>>")
else:
    print("No valid path found")