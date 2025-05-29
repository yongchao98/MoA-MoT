from collections import deque

# Define the initial state
start = 0
target = 9
initial_light = 'red'

# Define the operations
def press_A(number, light):
    return number * 2, 'green' if light == 'red' else 'red'

def press_B(number, light):
    if light == 'green':
        return number + 1, 'red'
    return number, light

def press_C(number, light):
    return number - 2, 'green' if light == 'red' else 'red'

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, initial_light, [])])
    visited = set()
    
    while queue:
        number, light, sequence = queue.popleft()
        
        if (number, light) in visited:
            continue
        visited.add((number, light))
        
        if number == target:
            return sequence
        
        # Try each button
        new_number, new_light = press_A(number, light)
        queue.append((new_number, new_light, sequence + ['A']))
        
        new_number, new_light = press_B(number, light)
        if new_number != number:  # Only add if B was successful
            queue.append((new_number, new_light, sequence + ['B']))
        
        new_number, new_light = press_C(number, light)
        queue.append((new_number, new_light, sequence + ['C']))

# Find the sequence
sequence = find_shortest_sequence(start, target)
print(" → ".join(sequence))