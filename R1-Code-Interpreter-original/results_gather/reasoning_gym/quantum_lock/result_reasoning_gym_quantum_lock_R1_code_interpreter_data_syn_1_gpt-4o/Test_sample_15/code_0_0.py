from collections import deque

# Define the initial state
initial_number = 0
initial_light = 'red'  # 'red' or 'green'
target_number = 24

# Define the button operations
def press_A(number, light):
    return number * 2, 'green' if light == 'red' else 'red'

def press_B(number, light):
    return number + 2, 'green' if light == 'red' else 'red'

def press_C(number, light):
    if light == 'green':
        return number * 2, 'red'
    return number, light

# BFS to find the shortest sequence
def find_shortest_sequence():
    queue = deque([(initial_number, initial_light, [])])  # (number, light, sequence)
    visited = set()
    
    while queue:
        number, light, sequence = queue.popleft()
        
        # Check if we reached the target
        if number == target_number:
            return sequence
        
        # Check if this state has been visited
        if (number, light) in visited:
            continue
        visited.add((number, light))
        
        # Try pressing each button
        new_number, new_light = press_A(number, light)
        queue.append((new_number, new_light, sequence + ['A']))
        
        new_number, new_light = press_B(number, light)
        queue.append((new_number, new_light, sequence + ['B']))
        
        new_number, new_light = press_C(number, light)
        if new_number != number:  # Only add if C was effective
            queue.append((new_number, new_light, sequence + ['C']))

# Find and print the shortest sequence
shortest_sequence = find_shortest_sequence()
print(' → '.join(shortest_sequence))