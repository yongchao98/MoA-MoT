def can_finish_courses(num_courses, prerequisites):
    from collections import defaultdict

    # Create an adjacency list for the graph
    graph = defaultdict(list)
    for course, prereq in prerequisites:
        graph[prereq].append(course)

    # States: 0 = unvisited, 1 = visiting, 2 = visited
    state = [0] * num_courses

    def has_cycle(v):
        if state[v] == 1:  # Found a cycle
            return True
        if state[v] == 2:  # Already visited
            return False

        # Mark the current node as visiting
        state[v] = 1
        for neighbor in graph[v]:
            if has_cycle(neighbor):
                return True

        # Mark the current node as visited
        state[v] = 2
        return False

    # Check each course for cycles
    for course in range(num_courses):
        if has_cycle(course):
            return False

    return True

# Number of courses
num_courses = 9
# List of prerequisites
prerequisites = [(5, 3), (2, 0), (8, 7), (8, 6), (7, 3), (1, 5), (6, 3), (3, 4), (2, 6), (0, 6)]

# Check if we can finish all courses
result = can_finish_courses(num_courses, prerequisites)
print(result)