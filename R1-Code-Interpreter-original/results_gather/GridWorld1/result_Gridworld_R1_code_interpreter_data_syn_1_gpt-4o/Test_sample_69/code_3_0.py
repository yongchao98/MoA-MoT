def dfs_find_path(current, path, visited_goals, goals, obstacles, adjacency):
    if visited_goals == goals:
        return path

    for neighbor in adjacency.get(current, []):
        if neighbor not in path and neighbor not in obstacles:
            new_path = path + [neighbor]
            new_visited_goals = visited_goals | {neighbor} if neighbor in goals else visited_goals
            result = dfs_find_path(neighbor, new_path, new_visited_goals, goals, obstacles, adjacency)
            if result:
                return result

    return None

# Find the path using DFS
path = dfs_find_path(initial_position, [initial_position], set(), goals, obstacles, adjacency)
print(path)