import json
from collections import deque, defaultdict
import heapq

def find_shortest_path(start, end, adjacency, obstacles):
    if start in obstacles or end in obstacles:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
            
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_all_pairs_shortest_paths(positions, adjacency, obstacles):
    paths = {}
    for i in positions:
        for j in positions:
            if i != j:
                path = find_shortest_path(i, j, adjacency, obstacles)
                if path:
                    paths[(i, j)] = path
    return paths

def solve_gridworld(start, goals, obstacles, adjacency):
    # Find all shortest paths between start and goals, and between all pairs of goals
    all_positions = [start] + goals
    all_paths = find_all_pairs_shortest_paths(all_positions, adjacency, obstacles)
    
    # Use a greedy approach to find a good order to visit goals
    current = start
    unvisited_goals = set(goals)
    final_path = [start]
    
    while unvisited_goals:
        # Find the closest unvisited goal
        best_dist = float('inf')
        best_goal = None
        best_path = None
        
        for goal in unvisited_goals:
            if (current, goal) in all_paths:
                path = all_paths[(current, goal)]
                if len(path) < best_dist:
                    best_dist = len(path)
                    best_goal = goal
                    best_path = path
        
        if best_goal is None:
            return None  # No valid path found
            
        # Add the path to this goal (excluding the start position)
        final_path.extend(best_path[1:])
        current = best_goal
        unvisited_goals.remove(best_goal)
    
    return final_path

# Input data
start = "C3,2"
goals = ['C2,1', 'C2,5', 'C1,5', 'C5,3', 'C5,1', 'C5,4']
obstacles = ['C5,2', 'C3,5', 'C2,4', 'C1,1', 'C5,5', 'C3,4', 'C3,1', 'C4,5']

# Define adjacency dictionary
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

# Solve and print result
path = solve_gridworld(start, goals, obstacles, adjacency)
print(f"<<<{json.dumps(path)}>>>")