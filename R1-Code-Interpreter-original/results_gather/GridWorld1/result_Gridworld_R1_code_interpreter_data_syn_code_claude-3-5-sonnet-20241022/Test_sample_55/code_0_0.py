import json
from collections import deque
import heapq

def heuristic(a, b):
    # Manhattan distance
    x1, y1 = map(int, a.replace('C', '').split(','))
    x2, y2 = map(int, b.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_path(came_from, current):
    path = [current]
    while current in came_from:
        current = came_from[current]
        path.append(current)
    return path[::-1]

def a_star(start, goal, adjacency, obstacles):
    frontier = [(0, start)]
    came_from = {}
    cost_so_far = {start: 0}
    
    while frontier:
        _, current = heapq.heappop(frontier)
        
        if current == goal:
            return get_path(came_from, current)
            
        for next_pos in adjacency[current]:
            if next_pos in obstacles:
                continue
                
            new_cost = cost_so_far[current] + 1
            if next_pos not in cost_so_far or new_cost < cost_so_far[next_pos]:
                cost_so_far[next_pos] = new_cost
                priority = new_cost + heuristic(next_pos, goal)
                heapq.heappush(frontier, (priority, next_pos))
                came_from[next_pos] = current
    
    return None

def find_complete_path(start, goals, adjacency, obstacles):
    current = start
    unvisited_goals = set(goals)
    path = [start]
    
    while unvisited_goals:
        # Find nearest unvisited goal
        min_dist = float('inf')
        nearest_goal = None
        best_path = None
        
        for goal in unvisited_goals:
            current_path = a_star(current, goal, adjacency, obstacles)
            if current_path is not None:
                dist = len(current_path)
                if dist < min_dist:
                    min_dist = dist
                    nearest_goal = goal
                    best_path = current_path
        
        if nearest_goal is None:
            return None  # No valid path found
            
        # Add path to nearest goal (excluding start position if not first goal)
        path.extend(best_path[1:])
        current = nearest_goal
        unvisited_goals.remove(nearest_goal)
    
    return path

# Main execution
start = "C6,6"
goals = ['C5,6', 'C4,4', 'C4,6', 'C2,5', 'C4,1', 'C1,3', 'C6,5']
obstacles = ['C3,2', 'C6,7', 'C1,1', 'C2,3', 'C5,7', 'C5,3', 'C1,6', 'C5,5', 'C3,3', 'C4,5']

# Parse adjacency from the given dictionary
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    # ... (full adjacency list as provided)
}

# Find the complete path
final_path = find_complete_path(start, goals, adjacency, obstacles)

# Output the result
if final_path:
    print(f"<<<{json.dumps(final_path)}>>>")
else:
    print("No valid path found")