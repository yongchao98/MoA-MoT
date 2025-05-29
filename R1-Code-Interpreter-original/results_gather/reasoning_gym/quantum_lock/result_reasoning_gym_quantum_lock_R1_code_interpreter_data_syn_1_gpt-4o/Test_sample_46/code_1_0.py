from collections import deque

# Define the initial state
initial_state = (0, 'red')  # (number, light color)
target = 14

# Define the operations
operations = {
    'A': lambda x: (x[0] + 2, 'green') if x[1] == 'red' else None,
    'B': lambda x: (x[0] * 2, 'red') if x[1] == 'green' else None,
    'C': lambda x: (x[0] + 1, 'green') if x[1] == 'red' else None
}

# BFS to find the shortest path
def bfs(initial_state, target):
    queue = deque([(initial_state, [])])  # (state, path)
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if we reached the target
        if current_state[0] == target:
            return path
        
        # Explore the next states
        for button, operation in operations.items():
            next_state = operation(current_state)
            if next_state and next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + [button]))

# Find the shortest sequence of button presses
shortest_path = bfs(initial_state, target)
print(" → ".join(shortest_path))