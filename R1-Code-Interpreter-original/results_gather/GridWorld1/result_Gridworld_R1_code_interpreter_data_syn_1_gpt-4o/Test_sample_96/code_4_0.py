import json
from heapq import heappop, heappush

# Define the grid and problem parameters
initial_position = "C2,6"
goals = {'C2,1', 'C2,3', 'C3,4', 'C5,2', 'C5,6', 'C4,3', 'C6,4'}
obstacles = {'C1,3', 'C5,5', 'C2,5', 'C5,3', 'C3,3', 'C6,5', 'C1,6', 'C1,2', 'C3,5', 'C6,2'}
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

# Heuristic function: Manhattan distance to the nearest goal
def heuristic(position, goals):
    x, y = map(int, position[1:].split(','))
    return min(abs(x - int(g[1])) + abs(y - int(g[3])) for g in goals)

# A* search to find a path visiting all goals
def a_star_path(initial, goals, obstacles, adjacency):
    open_set = []
    heappush(open_set, (0, initial, [initial], set()))
    while open_set:
        _, current, path, visited_goals = heappop(open_set)
        if visited_goals == goals:
            return path
        for neighbor in adjacency.get(current, []):
            if neighbor not in obstacles and neighbor not in path:
                new_visited_goals = visited_goals | {neighbor} if neighbor in goals else visited_goals
                cost = len(path) + heuristic(neighbor, goals - new_visited_goals)
                heappush(open_set, (cost, neighbor, path + [neighbor], new_visited_goals))
    return None

# Find the path using A* search
path = a_star_path(initial_position, goals, obstacles, adjacency)

# Output the path as JSON
print(json.dumps(path))