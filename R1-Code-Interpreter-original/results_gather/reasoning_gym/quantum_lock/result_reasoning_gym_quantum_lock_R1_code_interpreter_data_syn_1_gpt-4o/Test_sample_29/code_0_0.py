from collections import deque

# Define the initial state
start = (0, 'red')  # (number, light color)
target = 27

# Define the button operations
def press_A(number, light):
    if light == 'green':
        return number + 2, 'red'
    return number, light

def press_B(number, light):
    return number - 3, 'green' if light == 'red' else 'red'

def press_C(number, light):
    return number * 2, 'green' if light == 'red' else 'red'

# Breadth-first search to find the shortest sequence
queue = deque([(start, [])])
visited = set()

while queue:
    (current_number, current_light), path = queue.popleft()
    
    if current_number == target:
        print(" → ".join(path))
        break
    
    if (current_number, current_light) in visited:
        continue
    visited.add((current_number, current_light))
    
    # Try each button
    next_number, next_light = press_A(current_number, current_light)
    queue.append(((next_number, next_light), path + ['A']))
    
    next_number, next_light = press_B(current_number, current_light)
    queue.append(((next_number, next_light), path + ['B']))
    
    next_number, next_light = press_C(current_number, current_light)
    queue.append(((next_number, next_light), path + ['C']))