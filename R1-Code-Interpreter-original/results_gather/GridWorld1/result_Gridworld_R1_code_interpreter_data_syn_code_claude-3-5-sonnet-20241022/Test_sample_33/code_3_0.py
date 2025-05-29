from collections import deque
import json
import heapq

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_path(start, end, adjacency, obstacles):
    if start == end:
        return [start]
    
    frontier = [(manhattan_distance(start, end), start, [start])]
    visited = set()
    
    while frontier:
        _, current, path = heapq.heappop(frontier)
        
        if current == end:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                new_path = path + [next_pos]
                priority = len(new_path) + manhattan_distance(next_pos, end)
                heapq.heappush(frontier, (priority, next_pos, new_path))
    
    return None

def find_nearest_unvisited_goal(current, unvisited_goals, adjacency, obstacles, visited_goals):
    best_dist = float('inf')
    best_goal = None
    best_path = None
    
    for goal in unvisited_goals:
        path = get_path(current, goal, adjacency, obstacles)
        if path:
            dist = len(path) - 1  # Subtract 1 to get actual steps needed
            if dist < best_dist:
                best_dist = dist
                best_goal = goal
                best_path = path
    
    return best_goal, best_path

def find_complete_path(start, goals, adjacency, obstacles):
    current = start
    unvisited_goals = set(goals)
    visited_goals = set()
    final_path = [start]
    
    while unvisited_goals:
        next_goal, path = find_nearest_unvisited_goal(current, unvisited_goals, adjacency, obstacles, visited_goals)
        
        if not path:
            return None
            
        final_path.extend(path[1:])  # Skip the first position as it's already in the path
        current = next_goal
        unvisited_goals.remove(next_goal)
        visited_goals.add(next_goal)
    
    return final_path

# Problem setup
start = "C6,2"
goals = ['C6,3', 'C2,1', 'C3,6', 'C6,5', 'C2,2', 'C2,4', 'C6,4']
obstacles = ['C4,5', 'C1,2', 'C5,5', 'C3,5', 'C2,3', 'C6,6', 'C4,2', 'C5,4', 'C4,6', 'C5,6']

# Define adjacency dictionary
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"],
    "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5", "C5,5"],
    "C4,6": ["C4,5", "C3,6", "C5,6"],
    "C5,1": ["C5,2", "C4,1", "C6,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2", "C6,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3", "C6,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4", "C6,4"],
    "C5,5": ["C5,4", "C5,6", "C4,5", "C6,5"],
    "C5,6": ["C5,5", "C4,6", "C6,6"],
    "C6,1": ["C6,2", "C5,1"],
    "C6,2": ["C6,1", "C6,3", "C5,2"],
    "C6,3": ["C6,2", "C6,4", "C5,3"],
    "C6,4": ["C6,3", "C6,5", "C5,4"],
    "C6,5": ["C6,4", "C6,6", "C5,5"],
    "C6,6": ["C6,5", "C5,6"]
}

# Find solution path
solution = find_complete_path(start, goals, adjacency, obstacles)

# Format and print result
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No valid path found")