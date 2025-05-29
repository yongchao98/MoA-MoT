from collections import deque

# Define the initial state
start_number = 0
target_number = 43
start_light = 'red'

# Define the operations for each button
def apply_button_A(number):
    return number * 2

def apply_button_B(number, light):
    if light == 'red':
        return number - 2
    return None

def apply_button_C(number):
    return number + 3

# Define the BFS function
def bfs(start_number, target_number, start_light):
    # Queue for BFS: stores tuples of (current_number, current_light, path)
    queue = deque([(start_number, start_light, [])])
    # Set to keep track of visited states
    visited = set()
    
    while queue:
        current_number, current_light, path = queue.popleft()
        
        # Check if we have reached the target number
        if current_number == target_number:
            return path
        
        # If the current state has been visited, skip it
        if (current_number, current_light) in visited:
            continue
        
        # Mark the current state as visited
        visited.add((current_number, current_light))
        
        # Calculate the next states for each button
        # Button A
        next_number_A = apply_button_A(current_number)
        next_light_A = 'green' if current_light == 'red' else 'red'
        queue.append((next_number_A, next_light_A, path + ['A']))
        
        # Button B
        next_number_B = apply_button_B(current_number, current_light)
        if next_number_B is not None:
            next_light_B = 'green' if current_light == 'red' else 'red'
            queue.append((next_number_B, next_light_B, path + ['B']))
        
        # Button C
        next_number_C = apply_button_C(current_number)
        next_light_C = 'green' if current_light == 'red' else 'red'
        queue.append((next_number_C, next_light_C, path + ['C']))
    
    return None

# Find the shortest sequence of button presses
shortest_path = bfs(start_number, target_number, start_light)
print(' → '.join(shortest_path))