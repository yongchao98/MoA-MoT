def can_finish_courses(num_courses, prerequisites):
    from collections import defaultdict

    # Create a graph from the prerequisites
    graph = defaultdict(list)
    for course, prereq in prerequisites:
        graph[prereq].append(course)

    # Visited states: 0 = unvisited, 1 = visiting, 2 = visited
    visited = [0] * num_courses

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

# Given data
num_courses = 5
prerequisites = [(4, 0), (1, 0), (3, 0), (2, 0)]

# Check if we can finish all courses
result = can_finish_courses(num_courses, prerequisites)
print(result)