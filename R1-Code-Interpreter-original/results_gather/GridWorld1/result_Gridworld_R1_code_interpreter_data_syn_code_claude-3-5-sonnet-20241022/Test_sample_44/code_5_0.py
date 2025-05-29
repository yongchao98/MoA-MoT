import json
from collections import deque

def get_manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.split(',')[0][1:]), map(int, pos1.split(',')[1])
    x2, y2 = map(int, pos2.split(',')[0][1:]), map(int, pos2.split(',')[1])
    return abs(list(x1)[0] - list(x2)[0]) + abs(list(y1)[0] - list(y2)[0])

def find_direct_path(start, end, adjacency, obstacles, visited):
    queue = deque([(start, [start])])
    seen = {start}
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
            
        for next_pos in adjacency[current]:
            if (next_pos not in seen and 
                next_pos not in obstacles and 
                next_pos not in visited):
                seen.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_best_next_goal(current, unvisited_goals, adjacency, obstacles, visited):
    best_goal = None
    best_path = None
    best_score = float('inf')
    
    for goal in unvisited_goals:
        path = find_direct_path(current, goal, adjacency, obstacles, visited)
        if path:
            # Score based on path length and remaining goals accessibility
            score = len(path)
            remaining_goals = unvisited_goals - {goal}
            if remaining_goals:
                # Add heuristic for accessibility to other goals
                min_dist_to_others = min(get_manhattan_distance(goal, g) for g in remaining_goals)
                score += min_dist_to_others * 0.5
            
            if score < best_score:
                best_score = score
                best_path = path
                best_goal = goal
                
    return best_goal, best_path

def find_strict_path(start, goals, adjacency, obstacles):
    current = start
    unvisited_goals = set(goals)
    final_path = [start]
    visited = {start}
    
    while unvisited_goals:
        next_goal, path = find_best_next_goal(
            current,
            unvisited_goals,
            adjacency,
            obstacles,
            visited - {current}  # Allow current position to be used
        )
        
        if not path:
            return None
            
        # Add new positions to path and visited set
        final_path.extend(path[1:])
        visited.update(path[1:])
        current = next_goal
        unvisited_goals.remove(next_goal)
        
        # If we've reached all goals, stop
        if not unvisited_goals:
            break
    
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

solution = find_strict_path(start, goals, adjacency, obstacles)
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No valid path found")