def can_finish_courses(num_courses, prerequisites):
    from collections import defaultdict

    # Create an adjacency list for the graph
    graph = defaultdict(list)
    for course, prereq in prerequisites:
        graph[prereq].append(course)

    # States: 0 = unvisited, 1 = visiting, 2 = visited
    visited = [0] * num_courses

    def dfs(course):
        if visited[course] == 1:  # Found a cycle
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
num_courses = 6
# Prerequisites list
prerequisites = [(5, 1), (4, 5), (2, 5), (4, 0), (3, 0), (1, 4), (3, 5), (1, 5), (5, 0)]

# Check if we can finish all courses
result = can_finish_courses(num_courses, prerequisites)
print(result)