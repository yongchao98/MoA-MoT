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

# Number of courses
num_courses = 6
# Prerequisites list
prerequisites = [(0, 1), (3, 1), (5, 4), (5, 1), (4, 2), (2, 1), (2, 0)]

# Check if all courses can be finished
result = can_finish_courses(num_courses, prerequisites)
print(result)