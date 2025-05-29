from collections import deque

# Define the initial state
initial_state = (0, 'red')  # (number, light color)
target = 42

# Define the button operations
def press_A(number, light):
    return number - 3, 'green' if light == 'red' else 'red'

def press_B(number, light):
    return number + 2, 'green' if light == 'red' else 'red'

def press_C(number, light):
    if light == 'red':
        return number * 3, 'green'
    return number, light

# Breadth-first search to find the shortest sequence
queue = deque([(initial_state, [])])  # (state, sequence of buttons)
visited = set()

while queue:
    (current_number, current_light), sequence = queue.popleft()
    
    if current_number == target:
        print(' → '.join(sequence))
        break
    
    if (current_number, current_light) in visited:
        continue
    
    visited.add((current_number, current_light))
    
    # Try each button
    next_state_A = press_A(current_number, current_light)
    next_state_B = press_B(current_number, current_light)
    next_state_C = press_C(current_number, current_light)
    
    queue.append((next_state_A, sequence + ['A']))
    queue.append((next_state_B, sequence + ['B']))
    queue.append((next_state_C, sequence + ['C']))