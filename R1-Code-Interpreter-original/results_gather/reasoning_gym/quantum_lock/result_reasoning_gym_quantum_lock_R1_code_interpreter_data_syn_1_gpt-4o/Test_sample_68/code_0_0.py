from collections import deque

# Define the initial state
start = 0
target = 10
initial_light = 'red'

# Define the operations for each button
def press_A(number, light):
    return number * 2, 'green' if light == 'red' else 'red'

def press_B(number, light):
    if light == 'red':
        return number - 1, 'green'
    return number, light

def press_C(number, light):
    if light == 'red':
        return number + 1, 'green'
    return number, light

# Use a queue to perform a breadth-first search
queue = deque([(start, initial_light, [])])
visited = set()

while queue:
    number, light, path = queue.popleft()
    
    # Check if we reached the target
    if number == target:
        print(' → '.join(path))
        break
    
    # If this state has been visited, skip it
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