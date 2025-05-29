from collections import deque

# Define the initial state
start = (0, 'red')  # (number, light color)
target = 30

# Define the operations for each button
def press_A(number, light):
    if light == 'green':
        return number * 2, 'red'
    return None

def press_B(number, light):
    return number + 3, 'green' if light == 'red' else 'red'

def press_C(number, light):
    return number * 3, 'green' if light == 'red' else 'red'

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, [])])  # (state, path)
    visited = set()
    
    while queue:
        (current_number, current_light), path = queue.popleft()
        
        if current_number == target:
            return path
        
        if (current_number, current_light) in visited:
            continue
        
        visited.add((current_number, current_light))
        
        # Try each button
        if current_light == 'green':
            new_state = press_A(current_number, current_light)
            if new_state:
                queue.append((new_state, path + ['A']))
        
        new_state = press_B(current_number, current_light)
        queue.append((new_state, path + ['B']))
        
        new_state = press_C(current_number, current_light)
        queue.append((new_state, path + ['C']))
    
    return None

# Find the shortest sequence
sequence = find_shortest_sequence(start, target)
print(" → ".join(sequence))