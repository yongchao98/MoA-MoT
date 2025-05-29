from collections import deque

# Define the initial state
initial_state = (0, 'red')  # (number, light)

# Define the target number
target_number = 40

# Define the button operations
def press_button(state, button):
    number, light = state
    if button == 'A':
        if light == 'red':
            number += 3
        light = 'green' if light == 'red' else 'red'
    elif button == 'B':
        if light == 'green':
            number += 3
        light = 'green' if light == 'red' else 'red'
    elif button == 'C':
        number -= 2
        light = 'green' if light == 'red' else 'red'
    return (number, light)

# Use BFS to find the shortest sequence
def find_shortest_sequence(initial_state, target_number):
    queue = deque([(initial_state, [])])  # (state, sequence of buttons)
    visited = set()
    
    while queue:
        current_state, sequence = queue.popleft()
        
        # Check if we reached the target number
        if current_state[0] == target_number:
            return sequence
        
        # If already visited this state, skip it
        if current_state in visited:
            continue
        
        # Mark this state as visited
        visited.add(current_state)
        
        # Try all possible button presses
        for button in ['A', 'B', 'C']:
            new_state = press_button(current_state, button)
            new_sequence = sequence + [button]
            queue.append((new_state, new_sequence))

# Find the shortest sequence
shortest_sequence = find_shortest_sequence(initial_state, target_number)
print(" → ".join(shortest_sequence))