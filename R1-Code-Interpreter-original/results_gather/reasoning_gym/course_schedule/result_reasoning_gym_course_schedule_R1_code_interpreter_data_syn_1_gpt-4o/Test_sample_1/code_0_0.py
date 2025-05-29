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

# Given prerequisites
prerequisites = [(8, 0), (7, 5), (7, 1), (3, 5), (8, 4), (4, 2), (1, 0), (6, 7), (4, 0), (2, 8), (8, 1), (2, 3), (6, 5), (2, 5), (5, 4)]
num_courses = 9

# Check if all courses can be finished
result = can_finish_courses(num_courses, prerequisites)
print(result)