import json
import heapq

# Define the grid parameters
initial_position = "C6,6"
goals = {'C1,1', 'C3,4', 'C4,1', 'C4,4', 'C3,2', 'C1,4', 'C2,4'}
obstacles = {'C2,5', 'C5,4', 'C3,1', 'C6,1', 'C4,6', 'C6,4', 'C3,3', 'C2,2', 'C1,6', 'C3,6'}
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

# Helper function to calculate Manhattan distance
def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1[1:].split(','))
    x2, y2 = map(int, pos2[1:].split(','))
    return abs(x1 - x2) + abs(y1 - y2)

# A* search to find a path visiting all goals
def a_star_path(start, goals, obstacles, adjacency):
    open_set = []
    heapq.heappush(open_set, (0, start, [start], set()))
    while open_set:
        _, position, path, visited_goals = heapq.heappop(open_set)
        if visited_goals == goals:
            return path
        for neighbor in adjacency.get(position, []):
            if neighbor not in obstacles and neighbor not in path:
                new_visited_goals = visited_goals | {neighbor} if neighbor in goals else visited_goals
                heuristic = sum(manhattan_distance(neighbor, goal) for goal in goals - new_visited_goals)
                heapq.heappush(open_set, (len(path) + heuristic, neighbor, path + [neighbor], new_visited_goals))
    return None

# Find the path
path = a_star_path(initial_position, goals, obstacles, adjacency)

# Output the path as a JSON list
print(json.dumps(path))