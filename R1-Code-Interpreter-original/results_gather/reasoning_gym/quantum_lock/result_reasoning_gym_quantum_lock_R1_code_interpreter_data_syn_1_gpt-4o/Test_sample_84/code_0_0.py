from collections import deque

# Define the initial state
start_value = 0
target_value = 34
initial_light = 'red'

# Define the operations
def press_A(value, light):
    return value - 2, 'green' if light == 'red' else 'red'

def press_B(value, light):
    if light == 'red':
        return value + 2, 'green'
    return value, light

def press_C(value, light):
    return value * 3, 'green' if light == 'red' else 'red'

# BFS to find the shortest sequence
def find_shortest_sequence(start, target, initial_light):
    queue = deque([(start, initial_light, [])])  # (current_value, current_light, sequence)
    visited = set((start, initial_light))
    
    while queue:
        current_value, current_light, sequence = queue.popleft()
        
        if current_value == target:
            return sequence
        
        # Try each button
        # Button A
        new_value, new_light = press_A(current_value, current_light)
        if (new_value, new_light) not in visited:
            visited.add((new_value, new_light))
            queue.append((new_value, new_light, sequence + ['A']))
        
        # Button B
        new_value, new_light = press_B(current_value, current_light)
        if (new_value, new_light) not in visited:
            visited.add((new_value, new_light))
            queue.append((new_value, new_light, sequence + ['B']))
        
        # Button C
        new_value, new_light = press_C(current_value, current_light)
        if (new_value, new_light) not in visited:
            visited.add((new_value, new_light))
            queue.append((new_value, new_light, sequence + ['C']))
    
    return None

# Find the shortest sequence
sequence = find_shortest_sequence(start_value, target_value, initial_light)
print(" → ".join(sequence))