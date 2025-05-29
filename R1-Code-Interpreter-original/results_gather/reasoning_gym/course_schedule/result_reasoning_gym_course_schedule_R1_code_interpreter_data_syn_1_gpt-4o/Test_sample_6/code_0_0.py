def can_finish_courses(num_courses, prerequisites):
    from collections import defaultdict

    # Create an adjacency list for the graph
    graph = defaultdict(list)
    for course, prereq in prerequisites:
        graph[prereq].append(course)

    # Visited states: 0 = not visited, 1 = visiting, 2 = visited
    visited = [0] * num_courses

    def dfs(course):
        if visited[course] == 1:  # Cycle detected
            return False
        if visited[course] == 2:  # Already processed
            return True

        # Mark the course as visiting
        visited[course] = 1

        # Visit all the neighbors
        for neighbor in graph[course]:
            if not dfs(neighbor):
                return False

        # Mark the course as visited
        visited[course] = 2
        return True

    # Check each course
    for course in range(num_courses):
        if not dfs(course):
            return False

    return True

# Number of courses
num_courses = 10
# Prerequisites list
prerequisites = [(6, 2), (1, 3), (1, 4), (9, 0), (5, 7), (2, 5), (4, 7), (9, 1), (0, 3), (3, 7), (2, 4), (6, 3), (8, 1), (5, 4), (4, 8), (8, 4)]

# Check if all courses can be finished
result = can_finish_courses(num_courses, prerequisites)
print(result)