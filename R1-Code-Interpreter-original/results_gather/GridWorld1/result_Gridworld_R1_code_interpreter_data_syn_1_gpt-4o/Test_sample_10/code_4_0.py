import heapq

# Define the grid and problem parameters
initial_position = "C1,3"
goals = {'C5,1', 'C4,3', 'C5,4', 'C1,4', 'C3,4', 'C3,2'}
obstacles = {'C5,2', 'C5,3', 'C4,5', 'C5,5', 'C3,5', 'C2,4', 'C1,5', 'C2,5'}
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

# Heuristic function: Manhattan distance to the nearest unvisited goal
def heuristic(position, goals):
    x, y = map(int, position[1:].split(','))
    return min(abs(x - int(goal[1:goal.index(',')])) + abs(y - int(goal[goal.index(',')+1:])) for goal in goals)

# A* search to find a path visiting all goals
def a_star_search(initial, goals, obstacles, adjacency):
    open_set = []
    heapq.heappush(open_set, (0, initial, [initial], set()))
    visited = set()

    while open_set:
        _, current_position, path, visited_goals = heapq.heappop(open_set)

        # If all goals are visited, return the path
        if visited_goals == goals:
            return path

        # Explore neighbors
        for neighbor in adjacency.get(current_position, []):
            if neighbor not in obstacles and (neighbor, frozenset(visited_goals)) not in visited:
                new_visited_goals = visited_goals.copy()
                if neighbor in goals:
                    new_visited_goals.add(neighbor)
                new_path = path + [neighbor]
                cost = len(new_path) + heuristic(neighbor, goals - new_visited_goals)
                heapq.heappush(open_set, (cost, neighbor, new_path, new_visited_goals))
                visited.add((neighbor, frozenset(new_visited_goals)))

    return None

# Find the path
path = a_star_search(initial_position, goals, obstacles, adjacency)
print(path)