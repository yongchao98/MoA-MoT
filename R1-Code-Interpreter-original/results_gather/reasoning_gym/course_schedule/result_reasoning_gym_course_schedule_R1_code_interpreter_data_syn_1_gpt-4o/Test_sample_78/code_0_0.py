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
num_courses = 8
# Prerequisites list
prerequisites = [(0, 4), (1, 5), (0, 1), (6, 0), (6, 5), (4, 1), (0, 2), (1, 3), (3, 5), (3, 0), (2, 4), (7, 1)]

# Check if all courses can be finished
result = can_finish_courses(num_courses, prerequisites)
print(result)