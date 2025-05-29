from collections import deque

# Define the initial state
start = 0
target = 8
initial_light = 'red'

# Define the operations
def press_A(number):
    return number * 2

def press_B(number):
    return number + 1

def press_C(number):
    return number - 1

# Define the state as a tuple (current_number, light_color, sequence_of_buttons)
initial_state = (start, initial_light, [])

# Use a queue to perform a breadth-first search
queue = deque([initial_state])
visited = set()

while queue:
    current_number, light_color, sequence = queue.popleft()
    
    # Check if we reached the target
    if current_number == target:
        print(" → ".join(sequence))
        break
    
    # If this state has been visited, skip it
    if (current_number, light_color) in visited:
        continue
    visited.add((current_number, light_color))
    
    # Try pressing each button
    # Press A
    new_number = press_A(current_number)
    new_light = 'green' if light_color == 'red' else 'red'
    queue.append((new_number, new_light, sequence + ['A']))
    
    # Press B (only if light is green)
    if light_color == 'green':
        new_number = press_B(current_number)
        new_light = 'red'
        queue.append((new_number, new_light, sequence + ['B']))
    
    # Press C
    new_number = press_C(current_number)
    new_light = 'green' if light_color == 'red' else 'red'
    queue.append((new_number, new_light, sequence + ['C']))