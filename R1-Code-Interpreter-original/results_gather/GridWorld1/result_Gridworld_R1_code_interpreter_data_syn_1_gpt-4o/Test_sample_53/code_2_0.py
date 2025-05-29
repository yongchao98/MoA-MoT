def dfs_path(initial, goals, obstacles, adjacency):
    stack = [(initial, [initial], set())]
    while stack:
        current, path, visited_goals = stack.pop()
        if current in goals:
            visited_goals.add(current)
        if visited_goals == goals:
            return path
        for neighbor in adjacency.get(current, []):
            if neighbor not in obstacles and neighbor not in path:
                stack.append((neighbor, path + [neighbor], visited_goals.copy()))
    return []

# Find the path using DFS
path = dfs_path(initial_position, goals, obstacles, adjacency)
print(path)