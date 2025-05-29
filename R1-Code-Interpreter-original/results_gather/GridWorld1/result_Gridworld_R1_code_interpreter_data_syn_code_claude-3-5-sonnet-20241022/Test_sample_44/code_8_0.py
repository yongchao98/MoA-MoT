import json
from collections import deque

def find_path(start, target, adjacency, obstacles, must_visit=None):
    if must_visit is None:
        must_visit = set()
    
    queue = deque([(start, [start], set())])
    seen = {(start, frozenset())}
    
    while queue:
        current, path, visited_goals = queue.popleft()
        
        # If we reached the target and visited all required points
        if current == target and must_visit.issubset(visited_goals):
            return path
            
        for next_pos in adjacency[current]:
            if next_pos not in obstacles and next_pos not in path:
                new_visited = visited_goals.copy()
                if next_pos in must_visit:
                    new_visited.add(next_pos)
                
                state = (next_pos, frozenset(new_visited))
                if state not in seen:
                    seen.add(state)
                    queue.append((next_pos, path + [next_pos], new_visited))
    return None

def solve_gridworld(start, goals, adjacency, obstacles):
    # Try this specific sequence of goals
    sequence = [
        ("C4,4", "C3,4"),  # Move left first
        ("C3,4", "C3,3"),  # Move to C3,3 goal
        ("C3,3", "C3,2"),  # Move to C3,2 goal
        ("C3,2", "C4,2"),  # Move to C4,2 goal
        ("C4,2", "C4,1"),  # Move down
        ("C4,1", "C5,1"),  # Move to C5,1 goal
        ("C5,1", "C5,3"),  # Move to C5,3 goal
        ("C5,3", "C5,4")   # Finally reach C5,4 goal
    ]
    
    final_path = [start]
    current = start
    
    for start_pos, end_pos in sequence:
        if start_pos != current:
            continue
            
        # For each segment, find a path that might need to visit other goals along the way
        must_visit = set(g for g in goals if g not in final_path)
        path = find_path(start_pos, end_pos, adjacency, obstacles, must_visit)
        
        if path is None:
            return None
            
        final_path.extend(path[1:])  # Add all but the first position
        current = end_pos
        
        # Check if we've visited all goals
        if all(g in final_path for g in goals):
            return final_path
            
    return final_path if all(g in final_path for g in goals) else None

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

solution = solve_gridworld(start, goals, adjacency, obstacles)
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No valid path found")