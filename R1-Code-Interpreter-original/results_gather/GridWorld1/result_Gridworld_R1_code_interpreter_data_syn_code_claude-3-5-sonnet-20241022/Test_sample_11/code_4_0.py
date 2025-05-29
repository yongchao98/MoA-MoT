import json
from heapq import heappush, heappop
from typing import List, Dict, Set

def manhattan_distance(pos1: str, pos2: str) -> int:
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def find_path(start: str, goal: str, adjacency: Dict, obstacles: Set[str], must_visit: Set[str]) -> List[str]:
    if start == goal:
        return [start]
    
    queue = [(0, start, [start], set() & must_visit)]
    visited = {start: set() & must_visit}
    
    while queue:
        _, current, path, visited_goals = heappop(queue)
        
        if current == goal and must_visit.issubset(visited_goals):
            return path
            
        for next_pos in adjacency[current]:
            if next_pos not in obstacles:
                new_visited_goals = visited_goals | ({next_pos} & must_visit)
                if next_pos not in visited or not visited[next_pos].issubset(new_visited_goals):
                    visited[next_pos] = new_visited_goals
                    new_path = path + [next_pos]
                    priority = len(new_path) + manhattan_distance(next_pos, goal)
                    heappush(queue, (priority, next_pos, new_path, new_visited_goals))
    return []

def find_efficient_path(start: str, goals: List[str], adjacency: Dict, obstacles: Set[str]) -> List[str]:
    current = start
    remaining_goals = set(goals)
    path = [start]
    
    while remaining_goals:
        # Find the best next goal that minimizes total path length
        best_next_goal = None
        best_path_segment = None
        min_length = float('inf')
        
        for goal in remaining_goals:
            path_segment = find_path(current, goal, adjacency, obstacles, remaining_goals)
            if path_segment and len(path_segment) < min_length:
                min_length = len(path_segment)
                best_next_goal = goal
                best_path_segment = path_segment
        
        if best_next_goal is None:
            return []
        
        path.extend(best_path_segment[1:])
        current = best_next_goal
        remaining_goals.remove(best_next_goal)
    
    return path

start_pos = "C3,2"
goals = ['C1,4', 'C1,2', 'C3,5', 'C2,5', 'C2,1', 'C3,3', 'C2,2']
obstacles = {'C4,5', 'C4,6', 'C5,4', 'C2,3', 'C1,1', 'C5,6', 'C3,1', 'C5,3', 'C5,2', 'C4,3'}

adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"], "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"], "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"], "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"], "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"], "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"], "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1", "C5,1"], "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"], "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5", "C5,5"], "C4,6": ["C4,5", "C3,6", "C5,6"],
    "C5,1": ["C5,2", "C4,1"], "C5,2": ["C5,1", "C5,3", "C4,2"], "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"], "C5,5": ["C5,4", "C5,6", "C4,5"], "C5,6": ["C5,5", "C4,6"]
}

final_path = find_efficient_path(start_pos, goals, adjacency, obstacles)

if final_path and all(goal in final_path for goal in goals):
    print(f"<<<{json.dumps(final_path)}>>>")
else:
    print("No valid path found")