import json
from collections import deque
import heapq

def find_shortest_path(adjacency, start, end, obstacles):
    if start in obstacles or end in obstacles:
        return None
    
    queue = [(0, [start])]
    visited = {start}
    
    while queue:
        (dist, path) = heapq.heappop(queue)
        current = path[-1]
        
        if current == end:
            return path
            
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                visited.add(next_pos)
                new_path = list(path)
                new_path.append(next_pos)
                heapq.heappush(queue, (dist + 1, new_path))
    
    return None

def find_nearest_unvisited_goal(current, unvisited_goals, adjacency, obstacles):
    shortest_dist = float('inf')
    nearest_goal = None
    best_path = None
    
    for goal in unvisited_goals:
        path = find_shortest_path(adjacency, current, goal, obstacles)
        if path and len(path) < shortest_dist:
            shortest_dist = len(path)
            nearest_goal = goal
            best_path = path
    
    return nearest_goal, best_path

def find_complete_path(start, goals, obstacles, adjacency):
    current = start
    complete_path = [start]
    unvisited_goals = set(goals)
    
    while unvisited_goals:
        nearest_goal, path = find_nearest_unvisited_goal(current, unvisited_goals, adjacency, obstacles)
        if not nearest_goal:
            return None  # No valid path found
        
        # Add the path (excluding the start position which is already included)
        complete_path.extend(path[1:])
        current = nearest_goal
        unvisited_goals.remove(nearest_goal)
    
    return complete_path

# Define the problem parameters
start = "C4,2"
goals = ['C5,1', 'C4,6', 'C2,2', 'C3,1', 'C3,4', 'C2,4', 'C5,3']
obstacles = ['C3,5', 'C1,5', 'C2,1', 'C2,3', 'C4,1', 'C5,5', 'C2,6', 'C4,3', 'C1,1', 'C1,6']
adjacency = {
    # ... (your provided adjacency dict)
}

# Find the solution
solution = find_complete_path(start, goals, obstacles, adjacency)

# Format and print the solution
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No solution found")