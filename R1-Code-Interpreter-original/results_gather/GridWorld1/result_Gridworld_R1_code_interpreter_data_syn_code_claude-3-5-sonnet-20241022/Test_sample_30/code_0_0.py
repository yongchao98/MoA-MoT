import json
from collections import deque
import heapq

def bfs_path(graph, start, end, obstacles):
    if start in obstacles or end in obstacles:
        return None
    
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        current = path[-1]
        
        if current == end:
            return path
            
        for neighbor in graph[current]:
            if neighbor not in visited and neighbor not in obstacles:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    return None

def find_shortest_valid_path(graph, goals, start, obstacles):
    # Current state: position and unvisited goals
    initial_state = (start, frozenset(goals))
    
    # Priority queue for Dijkstra's algorithm
    pq = [(0, [start], initial_state)]
    visited = {initial_state}
    
    while pq:
        cost, path, (current, remaining_goals) = heapq.heappop(pq)
        
        # If we've visited all goals, we're done
        if not remaining_goals:
            return path
        
        # Try to visit each remaining goal
        for goal in remaining_goals:
            # Find path to this goal
            goal_path = bfs_path(graph, current, goal, obstacles)
            if goal_path is None:
                continue
                
            # Create new state
            new_path = path + goal_path[1:]  # Exclude start position to avoid duplicates
            new_remaining = remaining_goals - {goal}
            new_state = (goal, new_remaining)
            
            if new_state not in visited:
                visited.add(new_state)
                heapq.heappush(pq, (len(new_path), new_path, new_state))
    
    return None

# Define the problem
graph = {
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

start = "C2,5"
goals = ["C1,3", "C4,1", "C2,2", "C4,4", "C3,1", "C1,4"]
obstacles = ["C3,2", "C3,4", "C3,5", "C3,3", "C5,4", "C2,3", "C5,1", "C5,3"]

# Find solution
path = find_shortest_valid_path(graph, goals, start, obstacles)

# Output result
if path:
    print(f"<<<{json.dumps(path)}>>>")
else:
    print("No solution found")