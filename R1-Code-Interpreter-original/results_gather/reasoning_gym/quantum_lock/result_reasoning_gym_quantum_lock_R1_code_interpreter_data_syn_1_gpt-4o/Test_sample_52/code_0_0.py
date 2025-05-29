from collections import deque

# Define the initial state
start = 0
target = 18
initial_light = 'red'

# Define the button operations
def press_A(number, light):
    return number + 3, 'green' if light == 'red' else 'red'

def press_B(number, light):
    if light == 'green':
        return number * 2, 'red'
    return None

def press_C(number, light):
    if light == 'green':
        return number - 3, 'red'
    return None

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, initial_light, [])])  # (current number, light state, sequence of presses)
    visited = set()
    
    while queue:
        current, light, sequence = queue.popleft()
        
        if current == target:
            return sequence
        
        if (current, light) in visited:
            continue
        visited.add((current, light))
        
        # Try pressing each button
        new_number, new_light = press_A(current, light)
        queue.append((new_number, new_light, sequence + ['A']))
        
        if light == 'green':
            new_number, new_light = press_B(current, light)
            if new_number is not None:
                queue.append((new_number, new_light, sequence + ['B']))
            
            new_number, new_light = press_C(current, light)
            if new_number is not None:
                queue.append((new_number, new_light, sequence + ['C']))

# Find the shortest sequence
shortest_sequence = find_shortest_sequence(start, target)
print(' → '.join(shortest_sequence))