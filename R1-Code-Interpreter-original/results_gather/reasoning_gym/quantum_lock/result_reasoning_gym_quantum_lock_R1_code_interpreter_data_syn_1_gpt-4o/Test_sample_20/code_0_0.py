from collections import deque

# Define the initial state
start = 0
target = 19
initial_light = 'red'

# Define the operations
operations = {
    'A': lambda x, light: (x + 3, 'green' if light == 'red' else 'red') if light == 'green' else (x, light),
    'B': lambda x, light: (x - 3, 'green' if light == 'red' else 'red'),
    'C': lambda x, light: (x + 2, 'green' if light == 'red' else 'red') if light == 'red' else (x, light)
}

# Use BFS to find the shortest sequence
queue = deque([(start, initial_light, [])])
visited = set()

while queue:
    current_value, current_light, path = queue.popleft()
    
    if current_value == target:
        print(' → '.join(path))
        break
    
    for button, operation in operations.items():
        new_value, new_light = operation(current_value, current_light)
        if (new_value, new_light) not in visited:
            visited.add((new_value, new_light))
            queue.append((new_value, new_light, path + [button]))