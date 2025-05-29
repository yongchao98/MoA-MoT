from collections import deque

# Define the initial state
start = 0
target = 44
initial_light = 'red'

# Define the operations for each button
operations = {
    'A': lambda x, light: (x + 3, 'green' if light == 'red' else 'red'),
    'B': lambda x, light: (x + 3, 'green' if light == 'red' else 'red') if light == 'red' else (x, light),
    'C': lambda x, light: (x + 2, 'green' if light == 'red' else 'red') if light == 'green' else (x, light)
}

# BFS to find the shortest sequence
def find_shortest_sequence(start, target, initial_light):
    queue = deque([(start, initial_light, [])])
    visited = set()

    while queue:
        current_value, current_light, path = queue.popleft()

        if current_value == target:
            return path

        for button, operation in operations.items():
            new_value, new_light = operation(current_value, current_light)
            if (new_value, new_light) not in visited:
                visited.add((new_value, new_light))
                queue.append((new_value, new_light, path + [button]))

# Find the shortest sequence
sequence = find_shortest_sequence(start, target, initial_light)
print(" → ".join(sequence))