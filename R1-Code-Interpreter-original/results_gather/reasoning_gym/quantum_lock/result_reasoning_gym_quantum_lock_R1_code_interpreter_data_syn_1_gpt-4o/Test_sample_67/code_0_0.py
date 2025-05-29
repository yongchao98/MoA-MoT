from collections import deque

# Define the initial state
initial_state = (0, 'red')  # (number, light color)
target = 27

# Define the button operations
def press_button(state, button):
    number, light = state
    if button == 'A' and light == 'green':
        number -= 2
    elif button == 'B':
        number += 2
    elif button == 'C' and light == 'green':
        number += 3
    # Toggle the light
    light = 'green' if light == 'red' else 'red'
    return (number, light)

# BFS to find the shortest sequence
def find_shortest_sequence(initial_state, target):
    queue = deque([(initial_state, [])])  # (state, path)
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        
        if current_state[0] == target:
            return path
        
        if current_state in visited:
            continue
        
        visited.add(current_state)
        
        for button in ['A', 'B', 'C']:
            new_state = press_button(current_state, button)
            new_path = path + [button]
            queue.append((new_state, new_path))

# Find the shortest sequence
shortest_sequence = find_shortest_sequence(initial_state, target)
print(" → ".join(shortest_sequence))