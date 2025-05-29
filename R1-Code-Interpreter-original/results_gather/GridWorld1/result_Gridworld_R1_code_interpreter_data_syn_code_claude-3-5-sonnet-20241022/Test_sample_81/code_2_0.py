import json
from collections import deque

def find_path(start, end, adjacency, obstacles):
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
                new_path = list(path)
                new_path.append(next_pos)
                queue.append((next_pos, new_path))
    
    return None

def find_nearest_goal(current, goals, adjacency, obstacles, visited_goals):
    shortest_path = None
    shortest_length = float('inf')
    nearest_goal = None
    
    for goal in goals:
        if goal in visited_goals:
            continue
        path = find_path(current, goal, adjacency, obstacles)
        if path and len(path) < shortest_length:
            shortest_length = len(path)
            shortest_path = path
            nearest_goal = goal
            
    return nearest_goal, shortest_path

def solve_gridworld(start, goals, obstacles, adjacency):
    current = start
    final_path = [start]
    visited_goals = set()
    
    while len(visited_goals) < len(goals):
        next_goal, path = find_nearest_goal(current, goals, adjacency, obstacles, visited_goals)
        
        if not path:
            return None
        
        # Add the path (excluding the starting position)
        final_path.extend(path[1:])
        visited_goals.add(next_goal)
        current = next_goal
    
    return final_path

# Problem setup
start = "C4,2"
goals = ['C5,1', 'C4,6', 'C2,2', 'C3,1', 'C3,4', 'C2,4', 'C5,3']
obstacles = ['C3,5', 'C1,5', 'C2,1', 'C2,3', 'C4,1', 'C5,5', 'C2,6', 'C4,3', 'C1,1', 'C1,6']

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
  "C5,1": ["C5,2", "C4,1"],
  "C5,2": ["C5,1", "C5,3", "C4,2"],
  "C5,3": ["C5,2", "C5,4", "C4,3"],
  "C5,4": ["C5,3", "C5,5", "C4,4"],
  "C5,5": ["C5,4", "C5,6", "C4,5"],
  "C5,6": ["C5,5", "C4,6"]
}

# Execute solution
solution = solve_gridworld(start, goals, obstacles, adjacency)

# Output result
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No solution found")