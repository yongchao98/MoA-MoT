import json
from collections import deque
import heapq

def get_manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.split(',')[0][1:]), map(int, pos1.split(',')[1])
    x2, y2 = map(int, pos2.split(',')[0][1:]), map(int, pos2.split(',')[1])
    return abs(list(x1)[0] - list(x2)[0]) + abs(list(y1)[0] - list(y2)[0])

def find_path(adjacency, start, target, obstacles, visited=None):
    if visited is None:
        visited = set()
    
    queue = [(0, start, [start])]
    seen = {start}
    
    while queue:
        _, current, path = heapq.heappop(queue)
        
        if current == target:
            return path
        
        for next_pos in adjacency[current]:
            if next_pos not in seen and next_pos not in obstacles:
                seen.add(next_pos)
                new_path = path + [next_pos]
                priority = len(new_path) + get_manhattan_distance(next_pos, target)
                if next_pos in visited:
                    priority += 10  # Penalty for revisiting
                heapq.heappush(queue, (priority, next_pos, new_path))
    
    return None

def find_optimal_path(adjacency, start, goals, obstacles):
    remaining_goals = set(goals)
    current = start
    final_path = [start]
    visited = {start}
    
    while remaining_goals:
        best_next_goal = None
        best_path = None
        min_cost = float('inf')
        
        for goal in remaining_goals:
            path = find_path(adjacency, current, goal, obstacles, visited)
            if path:
                # Calculate cost considering path length and goal clustering
                base_cost = len(path)
                # Add heuristic for remaining goals
                remaining_cost = sum(get_manhattan_distance(goal, g) for g in remaining_goals if g != goal)
                total_cost = base_cost + remaining_cost/len(remaining_goals)
                
                if total_cost < min_cost:
                    min_cost = total_cost
                    best_path = path
                    best_next_goal = goal
        
        if not best_path:
            return None
        
        # Add path (excluding start position)
        for pos in best_path[1:]:
            final_path.append(pos)
            visited.add(pos)
        
        current = best_next_goal
        remaining_goals.remove(best_next_goal)
    
    return final_path

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

solution_path = find_optimal_path(adjacency, start, goals, obstacles)
if solution_path:
    print(f"<<<{json.dumps(solution_path)}>>>")
else:
    print("No valid path found")