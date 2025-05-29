from collections import deque

# Define the initial state
initial_state = (0, 'red')  # (number, light)

# Define the target number
target_number = 36

# Define the button operations
def press_A(number, light):
    return (number + 3, 'green' if light == 'red' else 'red')

def press_B(number, light):
    if light == 'green':
        return (number + 3, 'red')
    return None

def press_C(number, light):
    if light == 'red':
        return (number - 3, 'green')
    return None

# BFS to find the shortest sequence
def find_shortest_sequence():
    queue = deque([(initial_state, [])])  # (state, sequence of buttons)
    visited = set()
    
    while queue:
        (current_number, current_light), sequence = queue.popleft()
        
        # Check if we reached the target
        if current_number == target_number:
            return sequence
        
        # Mark the current state as visited
        visited.add((current_number, current_light))
        
        # Try pressing each button
        for button, operation in [('A', press_A), ('B', press_B), ('C', press_C)]:
            new_state = operation(current_number, current_light)
            if new_state and new_state not in visited:
                queue.append((new_state, sequence + [button]))

# Find and print the shortest sequence
shortest_sequence = find_shortest_sequence()
print(" → ".join(shortest_sequence))