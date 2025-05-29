def can_finish_courses(num_courses, prerequisites):
    from collections import defaultdict

    # Create an adjacency list for the graph
    graph = defaultdict(list)
    for dest, src in prerequisites:
        graph[src].append(dest)

    # Visited states: 0 = unvisited, 1 = visiting, 2 = visited
    visited = [0] * num_courses

    def dfs(course):
        if visited[course] == 1:  # cycle detected
            return False
        if visited[course] == 2:  # already visited
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
num_courses = 5
# Prerequisites list
prerequisites = [(2, 4), (3, 1), (2, 0), (1, 2), (3, 4), (4, 0), (2, 3), (0, 4), (1, 4), (0, 3)]

# Check if we can finish all courses
result = can_finish_courses(num_courses, prerequisites)
print(result)