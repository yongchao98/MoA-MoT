from collections import deque

# Define the initial state
start_number = 0
target_number = 31
initial_light = 'red'

# Define the button operations
def press_button(state, button):
    number, light = state
    if button == 'A' and light == 'red':
        number += 2
    elif button == 'B' and light == 'red':
        number -= 3
    elif button == 'C':
        number += 3
    # Toggle the light
    light = 'green' if light == 'red' else 'red'
    return number, light

# Perform BFS to find the shortest sequence of button presses
def find_shortest_sequence(start, target):
    queue = deque([(start, [])])  # (state, path)
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        current_number, current_light = current_state
        
        # Check if we've reached the target
        if current_number == target:
            return path
        
        # Check if we've already visited this state
        if current_state in visited:
            continue
        visited.add(current_state)
        
        # Try each button
        for button in ['A', 'B', 'C']:
            new_state = press_button(current_state, button)
            new_path = path + [button]
            queue.append((new_state, new_path))

# Find the shortest sequence
start_state = (start_number, initial_light)
shortest_sequence = find_shortest_sequence(start_state, target_number)
print(" → ".join(shortest_sequence))