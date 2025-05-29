from collections import deque

# Define the initial state
start = 0
target = 39
initial_light = 'red'

# Define the operations for each button
def press_A(number, light):
    if light == 'green':
        return number * 2, 'red'
    return number, light

def press_B(number, light):
    return number * 2, 'green' if light == 'red' else 'red'

def press_C(number, light):
    return number + 3, 'green' if light == 'red' else 'red'

# Use a queue to perform a breadth-first search
queue = deque([(start, initial_light, [])])
visited = set()

while queue:
    current_number, current_light, path = queue.popleft()
    
    if current_number == target:
        print(' → '.join(path))
        break
    
    if (current_number, current_light) in visited:
        continue
    
    visited.add((current_number, current_light))
    
    # Try each button
    new_number, new_light = press_A(current_number, current_light)
    queue.append((new_number, new_light, path + ['A']))
    
    new_number, new_light = press_B(current_number, current_light)
    queue.append((new_number, new_light, path + ['B']))
    
    new_number, new_light = press_C(current_number, current_light)
    queue.append((new_number, new_light, path + ['C']))