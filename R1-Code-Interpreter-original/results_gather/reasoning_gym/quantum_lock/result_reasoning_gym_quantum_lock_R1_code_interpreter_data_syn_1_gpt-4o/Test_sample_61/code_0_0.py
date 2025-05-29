from collections import deque

# Define the initial state
start = 0
target = 31
initial_light = 'red'

# Define the operations for each button
def press_A(number, light):
    if light == 'red':
        return number * 3, 'green'
    return number, light

def press_B(number, light):
    return number + 3, 'green' if light == 'red' else 'red'

def press_C(number, light):
    if light == 'green':
        return number - 2, 'red'
    return number, light

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, initial_light, [])])
    visited = set()
    
    while queue:
        number, light, path = queue.popleft()
        
        if number == target:
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
        queue.append((new_number, new_light, path + ['C']))
    
    return None

# Find and print the shortest sequence
sequence = find_shortest_sequence(start, target)
print(" → ".join(sequence))