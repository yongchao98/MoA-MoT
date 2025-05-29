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

        # Visit all the prerequisites
        for next_course in graph[course]:
            if not dfs(next_course):
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
num_courses = 7
# List of prerequisites
prerequisites = [(4, 6), (1, 3), (6, 3), (6, 2), (1, 5), (5, 2), (0, 5), (2, 3), (5, 6), (2, 0)]

# Check if all courses can be finished
result = can_finish_courses(num_courses, prerequisites)
print(result)