from collections import deque

# Define the initial state
start_number = 0
target_number = 31
initial_light = 'red'  # Light starts as red

# Define the button operations
def press_button_A(number, light):
    if light == 'green':
        return number - 2, 'red'
    return number, 'green'

def press_button_B(number, light):
    return number - 2, 'green' if light == 'red' else 'red'

def press_button_C(number, light):
    return number + 3, 'green' if light == 'red' else 'red'

# Initialize the queue with the starting state
queue = deque([(start_number, initial_light, [])])
visited = set((start_number, initial_light))

# Perform BFS
while queue:
    current_number, current_light, path = queue.popleft()
    
    # Check if we've reached the target
    if current_number == target_number:
        print(' → '.join(path))
        break
    
    # Explore the next states
    for button, operation in [('A', press_button_A), ('B', press_button_B), ('C', press_button_C)]:
        next_number, next_light = operation(current_number, current_light)
        if (next_number, next_light) not in visited:
            visited.add((next_number, next_light))
            queue.append((next_number, next_light, path + [button]))