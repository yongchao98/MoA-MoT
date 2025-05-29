def can_finish_courses(num_courses, prerequisites):
    from collections import defaultdict

    # Create an adjacency list for the graph
    graph = defaultdict(list)
    for course, prereq in prerequisites:
        graph[prereq].append(course)

    # Visited states: 0 = unvisited, 1 = visiting, 2 = visited
    visited = [0] * num_courses

    def dfs(course):
        if visited[course] == 1:  # Cycle detected
            return False
        if visited[course] == 2:  # Already visited
            return True

        # Mark the course as being visited
        visited[course] = 1

        # Visit all the neighbors
        for neighbor in graph[course]:
            if not dfs(neighbor):
                return False

        # Mark the course as fully visited
        visited[course] = 2
        return True

    # Check each course
    for course in range(num_courses):
        if not dfs(course):
            return False

    return True

# Number of courses
num_courses = 9
# List of prerequisites
prerequisites = [(7, 3), (1, 6), (7, 0), (8, 1), (4, 7), (2, 4), (1, 4), (3, 4), (5, 8), (6, 3), (4, 3), (5, 7), (4, 0), (3, 0)]

# Check if all courses can be finished
result = can_finish_courses(num_courses, prerequisites)
print(result)