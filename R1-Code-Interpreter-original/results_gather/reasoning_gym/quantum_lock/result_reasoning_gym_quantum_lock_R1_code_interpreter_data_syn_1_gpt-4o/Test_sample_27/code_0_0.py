from collections import deque

# Define the initial state
start = 0
target = 16
initial_light = 'red'

# Define the button operations
def press_A(number, light):
    return number + 2, 'green' if light == 'red' else 'red'

def press_B(number, light):
    if light == 'green':
        return number - 1, 'red'
    return number, light

def press_C(number, light):
    if light == 'green':
        return number * 2, 'red'
    return number, light

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, initial_light, [])])
    visited = set()
    
    while queue:
        number, light, path = queue.popleft()
        
        if (number, light) in visited:
            continue
        visited.add((number, light))
        
        if number == target:
            return path
        
        # Try each button
        new_number, new_light = press_A(number, light)
        queue.append((new_number, new_light, path + ['A']))
        
        new_number, new_light = press_B(number, light)
        if new_number != number:  # Only append if B was effective
            queue.append((new_number, new_light, path + ['B']))
        
        new_number, new_light = press_C(number, light)
        if new_number != number:  # Only append if C was effective
            queue.append((new_number, new_light, path + ['C']))

# Find the shortest sequence
sequence = find_shortest_sequence(start, target)
print(" → ".join(sequence))