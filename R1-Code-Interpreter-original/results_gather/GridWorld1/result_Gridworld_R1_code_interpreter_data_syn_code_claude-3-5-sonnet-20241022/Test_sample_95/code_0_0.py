from heapq import heappush, heappop
from collections import defaultdict
import json

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def find_path_to_all_goals(start, goals, obstacles, adjacency):
    def get_nearest_unvisited_goal(current, unvisited_goals):
        return min(unvisited_goals, key=lambda g: manhattan_distance(current, g))
    
    def a_star(start, target, visited_positions):
        frontier = [(0, start, [start])]
        came_from = {start: None}
        cost_so_far = {start: 0}
        
        while frontier:
            _, current, path = heappop(frontier)
            
            if current == target:
                return path
            
            for next_pos in adjacency[current]:
                if next_pos in obstacles:
                    continue
                    
                new_cost = cost_so_far[current] + 1
                
                if next_pos not in cost_so_far or new_cost < cost_so_far[next_pos]:
                    cost_so_far[next_pos] = new_cost
                    priority = new_cost + manhattan_distance(next_pos, target)
                    heappush(frontier, (priority, next_pos, path + [next_pos]))
                    came_from[next_pos] = current
        
        return None

    complete_path = [start]
    current = start
    remaining_goals = set(goals)
    visited_positions = {start}

    while remaining_goals:
        next_goal = get_nearest_unvisited_goal(current, remaining_goals)
        path = a_star(current, next_goal, visited_positions)
        
        if not path:
            return None
        
        # Add path (excluding the start position which is already included)
        complete_path.extend(path[1:])
        current = next_goal
        remaining_goals.remove(next_goal)
        visited_positions.update(path)

    return complete_path

# Problem setup
start = "C2,2"
goals = ['C2,4', 'C3,2', 'C5,3', 'C2,3', 'C5,6', 'C1,4', 'C5,5']
obstacles = ['C1,2', 'C1,3', 'C4,1', 'C5,2', 'C3,1', 'C2,6', 'C4,5', 'C4,2', 'C4,6', 'C1,6']

# Load adjacency data
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

# Find path
path = find_path_to_all_goals(start, goals, obstacles, adjacency)

# Output result
print(f"<<<{json.dumps(path)}>>>")