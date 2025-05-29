from collections import deque

# Initial state
start_number = 0
target_number = 14
initial_light = 'red'  # Light starts as red

# Define the button operations
def press_A(number, light):
    if light == 'green':
        return number + 2, 'red'
    return number, light

def press_B(number, light):
    return number * 2, 'green' if light == 'red' else 'red'

def press_C(number, light):
    if light == 'green':
        return number + 1, 'red'
    return number, light

# BFS to find the shortest sequence
queue = deque([(start_number, initial_light, [])])
visited = set()

while queue:
    current_number, current_light, sequence = queue.popleft()
    
    # Check if we reached the target
    if current_number == target_number:
        print(' → '.join(sequence))
        break
    
    # Check if this state has been visited
    if (current_number, current_light) in visited:
        continue
    visited.add((current_number, current_light))
    
    # Try each button
    new_number, new_light = press_A(current_number, current_light)
    queue.append((new_number, new_light, sequence + ['A']))
    
    new_number, new_light = press_B(current_number, current_light)
    queue.append((new_number, new_light, sequence + ['B']))
    
    new_number, new_light = press_C(current_number, current_light)
    queue.append((new_number, new_light, sequence + ['C']))