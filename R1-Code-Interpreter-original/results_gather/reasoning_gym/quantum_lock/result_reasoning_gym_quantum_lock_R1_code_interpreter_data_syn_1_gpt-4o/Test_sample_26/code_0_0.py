from collections import deque

# Define the initial state
start = 0
target = 38
initial_light = 'red'

# Define the button operations
def press_button(state, button, light):
    number, light = state
    if button == 'A' and light == 'green':
        number *= 3
    elif button == 'B' and light == 'green':
        number += 2
    elif button == 'C':
        number += 3
    # Toggle the light
    light = 'green' if light == 'red' else 'red'
    return number, light

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([((start, initial_light), [])])
    visited = set()
    
    while queue:
        (current_number, current_light), path = queue.popleft()
        
        if current_number == target:
            return path
        
        if (current_number, current_light) in visited:
            continue
        
        visited.add((current_number, current_light))
        
        for button in ['A', 'B', 'C']:
            new_state = press_button((current_number, current_light), button, current_light)
            new_path = path + [button]
            queue.append((new_state, new_path))
    
    return None

# Find the shortest sequence
sequence = find_shortest_sequence(start, target)
print(" → ".join(sequence))