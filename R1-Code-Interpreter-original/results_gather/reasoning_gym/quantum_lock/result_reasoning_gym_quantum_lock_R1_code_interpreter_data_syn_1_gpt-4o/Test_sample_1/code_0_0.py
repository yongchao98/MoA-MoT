from collections import deque

# Initial state
start_number = 0
target_number = 46
initial_light = 'red'  # Light starts as red

# Define the operations for each button
def press_A(number, light):
    return number * 3, 'green' if light == 'red' else 'red'

def press_B(number, light):
    if light == 'green':
        return number - 2, 'red'
    return None

def press_C(number, light):
    if light == 'red':
        return number + 3, 'green'
    return None

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, initial_light, [])])  # (current number, light state, sequence of presses)
    visited = set()
    
    while queue:
        current_number, current_light, sequence = queue.popleft()
        
        # Check if we reached the target
        if current_number == target:
            return sequence
        
        # Avoid revisiting the same state
        if (current_number, current_light) in visited:
            continue
        visited.add((current_number, current_light))
        
        # Try each button
        # Button A
        new_number, new_light = press_A(current_number, current_light)
        queue.append((new_number, new_light, sequence + ['A']))
        
        # Button B
        result = press_B(current_number, current_light)
        if result:
            new_number, new_light = result
            queue.append((new_number, new_light, sequence + ['B']))
        
        # Button C
        result = press_C(current_number, current_light)
        if result:
            new_number, new_light = result
            queue.append((new_number, new_light, sequence + ['C']))
    
    return None

# Find the shortest sequence
sequence = find_shortest_sequence(start_number, target_number)
print(" → ".join(sequence))