import json
from collections import deque, defaultdict
import heapq

def find_shortest_path(adjacency, start, end, obstacles):
    if start in obstacles or end in obstacles:
        return None
    
    queue = [(0, start, [start])]
    visited = set()
    
    while queue:
        (cost, current, path) = heapq.heappop(queue)
        
        if current == end:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                new_path = path + [next_pos]
                heapq.heappush(queue, (cost + 1, next_pos, new_path))
    
    return None

def find_complete_path(adjacency, start, goals, obstacles):
    current_pos = start
    complete_path = [start]
    unvisited_goals = set(goals)
    visited = set([start])
    
    # If start is a goal, mark it as visited
    if start in unvisited_goals:
        unvisited_goals.remove(start)
    
    while unvisited_goals:
        # Find the nearest unvisited goal
        best_path = None
        best_goal = None
        min_length = float('inf')
        
        for goal in unvisited_goals:
            path = find_shortest_path(adjacency, current_pos, goal, obstacles)
            if path and len(path) < min_length:
                min_length = len(path)
                best_path = path
                best_goal = goal
        
        if not best_path:
            return None
        
        # Add the path to the goal (excluding the start position)
        for pos in best_path[1:]:
            if pos not in visited or pos in unvisited_goals:
                complete_path.append(pos)
                visited.add(pos)
        
        unvisited_goals.remove(best_goal)
        current_pos = best_goal
    
    return complete_path

# Define the problem parameters
start = "C3,1"
goals = ['C1,5', 'C3,4', 'C3,3', 'C4,4', 'C2,1', 'C2,5', 'C4,3']
obstacles = ['C3,2', 'C5,6', 'C5,4', 'C2,3', 'C3,5', 'C1,6', 'C5,1', 'C2,4', 'C1,1', 'C5,2']

# Define adjacency
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
  "C5,1": ["C5,2", "C4,1"],
  "C5,2": ["C5,1", "C5,3", "C4,2"],
  "C5,3": ["C5,2", "C5,4", "C4,3"],
  "C5,4": ["C5,3", "C5,5", "C4,4"],
  "C5,5": ["C5,4", "C5,6", "C4,5"],
  "C5,6": ["C5,5", "C4,6"]
}

# Find the complete path
path = find_complete_path(adjacency, start, goals, obstacles)

# Format the output
if path:
    output = json.dumps(path)
    print(f"<<<{output}>>>")
else:
    print("<<<No valid path found>>>")