from collections import deque

# Define the initial state
initial_state = (0, 'red')  # (number, light color)
target = 21

# Define the button operations
def press_A(number, light):
    if light == 'green':
        return number - 2, 'red'
    return number, light

def press_B(number, light):
    return number * 2, 'green' if light == 'red' else 'red'

def press_C(number, light):
    if light == 'green':
        return number + 3, 'red'
    return number, light

# Define the BFS function to find the shortest sequence
def find_shortest_sequence(initial_state, target):
    queue = deque([(initial_state, [])])  # (state, sequence of buttons)
    visited = set()
    
    while queue:
        (current_number, current_light), sequence = queue.popleft()
        
        if current_number == target:
            return sequence
        
        if (current_number, current_light) in visited:
            continue
        
        visited.add((current_number, current_light))
        
        # Try each button
        new_state_A = press_A(current_number, current_light)
        new_state_B = press_B(current_number, current_light)
        new_state_C = press_C(current_number, current_light)
        
        queue.append((new_state_A, sequence + ['A']))
        queue.append((new_state_B, sequence + ['B']))
        queue.append((new_state_C, sequence + ['C']))
    
    return None

# Find the shortest sequence
sequence = find_shortest_sequence(initial_state, target)
print(" → ".join(sequence))