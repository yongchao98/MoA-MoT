from collections import deque

# Define the initial state
start_state = (0, 'red')  # (number, light color)
target = 27

# Define the button operations
def apply_button(state, button):
    number, light = state
    if button == 'A' and light == 'green':
        number -= 3
    elif button == 'B' and light == 'green':
        number += 3
    elif button == 'C':
        number *= 3
    # Toggle the light
    light = 'green' if light == 'red' else 'red'
    return (number, light)

# Perform BFS to find the shortest sequence
def find_shortest_sequence(start_state, target):
    queue = deque([(start_state, [])])  # (state, path)
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if we've reached the target
        if current_state[0] == target:
            return path
        
        # Mark the current state as visited
        visited.add(current_state)
        
        # Try each button
        for button in ['A', 'B', 'C']:
            new_state = apply_button(current_state, button)
            if new_state not in visited:
                queue.append((new_state, path + [button]))

# Find the shortest sequence
shortest_sequence = find_shortest_sequence(start_state, target)
print(" → ".join(shortest_sequence))