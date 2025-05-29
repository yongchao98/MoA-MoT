import json
from heapq import heappush, heappop
from collections import defaultdict

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def find_nearest_unvisited_goal(current, goals, visited_goals):
    unvisited = [g for g in goals if g not in visited_goals]
    if not unvisited:
        return None
    return min(unvisited, key=lambda g: manhattan_distance(current, g))

def find_path(start, goals, obstacles, adjacency):
    goals = set(goals)
    obstacles = set(obstacles)
    
    def heuristic(pos, unvisited_goals):
        if not unvisited_goals:
            return 0
        return min(manhattan_distance(pos, g) for g in unvisited_goals)
    
    # Priority queue entries: (priority, path_length, current_pos, path, visited_goals)
    start_state = (0, 0, start, [start], frozenset())
    frontier = [start_state]
    seen = set()  # (current_pos, visited_goals)
    
    while frontier:
        _, path_len, current, path, visited_goals = heappop(frontier)
        
        state = (current, visited_goals)
        if state in seen:
            continue
        seen.add(state)
        
        # Update visited goals
        new_visited = set(visited_goals)
        if current in goals:
            new_visited.add(current)
            
        # Check if all goals are visited
        if len(new_visited) == len(goals):
            return path
            
        # Try all possible moves
        for next_pos in adjacency[current]:
            if next_pos in obstacles:
                continue
                
            new_path = path + [next_pos]
            unvisited = goals - new_visited
            
            # Calculate priority based on path length and heuristic
            h = heuristic(next_pos, unvisited)
            priority = path_len + 1 + h
            
            new_state = (priority, path_len + 1, next_pos, new_path, frozenset(new_visited))
            heappush(frontier, new_state)
    
    return None

# Problem setup
start = "C3,3"
goals = ['C2,6', 'C4,4', 'C1,3', 'C5,4', 'C2,3', 'C3,4', 'C5,5']
obstacles = ['C1,1', 'C6,4', 'C5,6', 'C1,4', 'C6,2', 'C4,3', 'C1,2', 'C4,1', 'C4,2', 'C4,6']

# Load adjacency from the provided dictionary (assume it's defined as 'adjacency')
adjacency = {
    # ... (your provided adjacency dictionary)
}

# Find path
path = find_path(start, goals, obstacles, adjacency)

# Output result
if path:
    result = json.dumps(path)
    print(f"<<<{result}>>>")
else:
    print("No valid path found")