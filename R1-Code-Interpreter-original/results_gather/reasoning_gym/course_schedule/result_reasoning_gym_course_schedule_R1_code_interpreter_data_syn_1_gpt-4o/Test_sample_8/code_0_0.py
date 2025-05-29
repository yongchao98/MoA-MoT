def can_finish_courses(num_courses, prerequisites):
    from collections import defaultdict

    # Create an adjacency list for the graph
    graph = defaultdict(list)
    for course, prereq in prerequisites:
        graph[prereq].append(course)

    # Visited set to keep track of visited nodes
    visited = set()
    # Recursion stack to detect cycles
    rec_stack = set()

    def dfs(course):
        if course in rec_stack:
            return False  # Cycle detected
        if course in visited:
            return True  # Already visited and no cycle found

        # Mark the course as visited and add to recursion stack
        visited.add(course)
        rec_stack.add(course)

        # Visit all the neighbors
        for neighbor in graph[course]:
            if not dfs(neighbor):
                return False

        # Remove from recursion stack
        rec_stack.remove(course)
        return True

    # Check each course
    for course in range(num_courses):
        if course not in visited:
            if not dfs(course):
                return False

    return True

# Number of courses
num_courses = 8
# List of prerequisites
prerequisites = [(5, 7), (0, 1), (2, 5), (6, 7), (3, 7), (0, 7), (2, 6), (4, 7), (1, 7), (5, 0), (3, 5), (4, 5), (1, 2), (6, 4)]

# Check if all courses can be finished
result = can_finish_courses(num_courses, prerequisites)
print(result)