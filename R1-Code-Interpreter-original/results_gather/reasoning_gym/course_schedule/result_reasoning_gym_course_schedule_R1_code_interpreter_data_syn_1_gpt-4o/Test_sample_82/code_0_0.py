def can_finish_courses(num_courses, prerequisites):
    from collections import defaultdict

    # Create a graph from the prerequisites
    graph = defaultdict(list)
    for course, prereq in prerequisites:
        graph[prereq].append(course)

    # Track visited nodes to detect cycles
    visited = [0] * num_courses  # 0 = unvisited, 1 = visiting, 2 = visited

    def dfs(course):
        if visited[course] == 1:  # Cycle detected
            return False
        if visited[course] == 2:  # Already visited
            return True

        # Mark the course as being visited
        visited[course] = 1
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
# List of prerequisites
prerequisites = [(0, 1), (5, 1), (4, 5), (3, 4), (2, 1), (4, 2), (2, 5), (0, 3), (5, 0), (3, 5)]

# Check if all courses can be finished
result = can_finish_courses(num_courses, prerequisites)
print(result)