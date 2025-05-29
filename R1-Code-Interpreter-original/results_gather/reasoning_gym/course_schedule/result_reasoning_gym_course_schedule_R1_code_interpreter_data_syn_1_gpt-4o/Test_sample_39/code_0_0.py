from collections import defaultdict, deque

def can_finish_courses(num_courses, prerequisites):
    # Create an adjacency list for the graph
    graph = defaultdict(list)
    indegree = [0] * num_courses
    
    # Build the graph and calculate indegrees
    for dest, src in prerequisites:
        graph[src].append(dest)
        indegree[dest] += 1
    
    # Initialize a queue with all nodes having indegree 0
    queue = deque([i for i in range(num_courses) if indegree[i] == 0])
    
    # Number of courses that can be completed
    completed_courses = 0
    
    while queue:
        course = queue.popleft()
        completed_courses += 1
        
        # Decrease the indegree of neighboring nodes
        for neighbor in graph[course]:
            indegree[neighbor] -= 1
            # If indegree becomes 0, add to queue
            if indegree[neighbor] == 0:
                queue.append(neighbor)
    
    # If completed courses equal total courses, return True
    return completed_courses == num_courses

# Given number of courses and prerequisites
num_courses = 7
prerequisites = [(0, 1), (2, 1), (3, 1), (5, 3), (4, 0), (5, 1), (2, 3), (6, 0)]

# Check if all courses can be finished
result = can_finish_courses(num_courses, prerequisites)
print(result)