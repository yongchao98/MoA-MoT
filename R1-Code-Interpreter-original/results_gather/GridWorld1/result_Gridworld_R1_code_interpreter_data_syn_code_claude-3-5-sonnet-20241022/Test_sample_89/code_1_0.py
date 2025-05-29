import json
from heapq import heappush, heappop

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

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

# Complete adjacency dictionary
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

# Find path
path = find_path(start, goals, obstacles, adjacency)

# Output result
if path:
    result = json.dumps(path)
    print(f"<<<{result}>>>")
else:
    print("No valid path found")