from collections import deque

# Define the initial state
initial_state = (0, 'red')  # (number, light color)
target = 15

# Define the button operations
def press_A(number, light):
    return number - 2, 'green' if light == 'red' else 'red'

def press_B(number, light):
    return number + 1, 'green' if light == 'red' else 'red'

def press_C(number, light):
    if light == 'green':
        return number * 2, 'red'
    return None  # C cannot be pressed when the light is red

# BFS to find the shortest sequence
def find_shortest_sequence(initial_state, target):
    queue = deque([(initial_state, [])])  # (state, sequence of buttons)
    visited = set()
    
    while queue:
        (current_number, current_light), sequence = queue.popleft()
        
        # Check if we reached the target
        if current_number == target:
            return sequence
        
        # If already visited this state, skip it
        if (current_number, current_light) in visited:
            continue
        visited.add((current_number, current_light))
        
        # Try each button
        # Press A
        new_state = press_A(current_number, current_light)
        queue.append((new_state, sequence + ['A']))
        
        # Press B
        new_state = press_B(current_number, current_light)
        queue.append((new_state, sequence + ['B']))
        
        # Press C (only if applicable)
        if current_light == 'green':
            new_state = press_C(current_number, current_light)
            if new_state:
                queue.append((new_state, sequence + ['C']))

# Find the shortest sequence
shortest_sequence = find_shortest_sequence(initial_state, target)
print(" → ".join(shortest_sequence))