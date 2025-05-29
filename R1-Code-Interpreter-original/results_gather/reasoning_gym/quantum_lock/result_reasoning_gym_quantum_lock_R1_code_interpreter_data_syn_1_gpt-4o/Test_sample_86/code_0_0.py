from collections import deque

# Define the initial state
start = 0
target = 12
initial_light = 'red'

# Define the operations for each button
operations = {
    'A': lambda x, light: (x - 1, 'red' if light == 'green' else 'green') if light == 'green' else (x, light),
    'B': lambda x, light: (x - 1, 'red' if light == 'green' else 'green'),
    'C': lambda x, light: (x + 2, 'red' if light == 'green' else 'green')
}

# Use a queue to perform a breadth-first search
queue = deque([(start, initial_light, [])])
visited = set()

while queue:
    current_value, current_light, path = queue.popleft()
    
    # Check if we have reached the target
    if current_value == target:
        print(' → '.join(path))
        break
    
    # Try each button
    for button, operation in operations.items():
        new_value, new_light = operation(current_value, current_light)
        
        # Avoid revisiting the same state
        if (new_value, new_light) not in visited:
            visited.add((new_value, new_light))
            queue.append((new_value, new_light, path + [button]))