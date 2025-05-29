from collections import deque

# Initial state
start_number = 0
start_light = 'red'
target_number = 24

# Define the button operations
def press_A(number, light):
    return number * 2, 'green' if light == 'red' else 'red'

def press_B(number, light):
    return number + 2, 'green' if light == 'red' else 'red'

def press_C(number, light):
    if light == 'green':
        return number + 2, 'red'
    return number, light

# BFS to find the shortest sequence
def find_shortest_sequence(start_number, start_light, target_number):
    queue = deque([(start_number, start_light, [])])
    visited = set()
    
    while queue:
        number, light, path = queue.popleft()
        
        if number == target_number:
            return path
        
        if (number, light) in visited:
            continue
        
        visited.add((number, light))
        
        # Try each button
        new_number, new_light = press_A(number, light)
        queue.append((new_number, new_light, path + ['A']))
        
        new_number, new_light = press_B(number, light)
        queue.append((new_number, new_light, path + ['B']))
        
        new_number, new_light = press_C(number, light)
        if new_number != number:  # Only add if C was effective
            queue.append((new_number, new_light, path + ['C']))
    
    return None

# Find the shortest sequence
sequence = find_shortest_sequence(start_number, start_light, target_number)
print(" → ".join(sequence))