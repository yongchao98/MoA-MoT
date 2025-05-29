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
    
    frontier = [(manhattan_distance(start, end), 0, start, [start])]
    visited = set()
    
    while frontier:
        _, cost, current, path = heapq.heappop(frontier)
        
        if current == end:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                new_path = path + [next_pos]
                new_cost = cost + 1
                priority = new_cost + manhattan_distance(next_pos, end) * 2
                heapq.heappush(frontier, (priority, new_cost, next_pos, new_path))
    
    return None

def find_complete_path(start, goals, adjacency, obstacles):
    current = start
    unvisited_goals = set(goals)
    path = [start]
    visited_positions = {start}
    
    while unvisited_goals:
        # Find best next goal considering both distance and path efficiency
        best_score = float('inf')
        best_goal = None
        best_path = None
        
        for goal in unvisited_goals:
            temp_path = get_path(current, goal, adjacency, obstacles)
            if temp_path:
                # Score based on path length and how many unvisited goals are near the path
                score = len(temp_path)
                best_score = score
                best_goal = goal
                best_path = temp_path
                break
        
        if best_path is None:
            return None
        
        # Add path segment (excluding start position which is already in path)
        path.extend(best_path[1:])
        current = best_goal
        unvisited_goals.remove(best_goal)
        
    return path

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