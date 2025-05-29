def can_finish_courses(num_courses, prerequisites):
    from collections import defaultdict

    # Create a graph from the prerequisites
    graph = defaultdict(list)
    for course, prereq in prerequisites:
        graph[prereq].append(course)

    # Track visited nodes to detect cycles
    visited = [0] * num_courses  # 0 = unvisited, 1 = visiting, 2 = visited

    def has_cycle(v):
        if visited[v] == 1:  # Found a cycle
            return True
        if visited[v] == 2:  # Already visited node
            return False

        # Mark the node as visiting
        visited[v] = 1
        for neighbor in graph[v]:
            if has_cycle(neighbor):
                return True

        # Mark the node as visited
        visited[v] = 2
        return False

    # Check each course for cycles
    for course in range(num_courses):
        if has_cycle(course):
            return False

    return True

# Number of courses
num_courses = 5
# List of prerequisites
prerequisites = [(1, 3), (0, 2), (1, 0), (3, 4), (3, 0), (4, 0), (2, 4), (2, 1)]

# Check if all courses can be finished
result = can_finish_courses(num_courses, prerequisites)
print(result)