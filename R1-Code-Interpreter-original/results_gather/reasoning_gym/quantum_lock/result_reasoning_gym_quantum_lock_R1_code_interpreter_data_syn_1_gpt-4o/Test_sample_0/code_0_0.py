from collections import deque

# Define the initial state
start = 0
target = 9

# Define the operations for each button
operations = {
    'A': lambda x: x + 1,
    'B': lambda x: x + 2,
    'C': lambda x: x - 1
}

# Use a queue to perform BFS
queue = deque([(start, [])])

# Perform BFS
while queue:
    current, path = queue.popleft()
    
    # Check if we've reached the target
    if current == target:
        print(" → ".join(path))
        break
    
    # Explore each button press
    for button, operation in operations.items():
        new_value = operation(current)
        new_path = path + [button]
        queue.append((new_value, new_path))