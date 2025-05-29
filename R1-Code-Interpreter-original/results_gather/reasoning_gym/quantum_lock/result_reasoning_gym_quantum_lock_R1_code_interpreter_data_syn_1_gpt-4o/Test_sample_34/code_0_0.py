from collections import deque

# Define the initial state
initial_state = (0, 'red')  # (number, light)

# Define the target number
target = 8

# Define the operations for each button
def press_A(number, light):
    return (number * 2, 'green' if light == 'red' else 'red')

def press_B(number, light):
    return (number - 1, 'green' if light == 'red' else 'red')

def press_C(number, light):
    if light == 'green':
        return (number + 1, 'red')
    return None

# Breadth-first search to find the shortest sequence
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
        
        # Try pressing each button
        next_state_A = press_A(current_number, current_light)
        next_state_B = press_B(current_number, current_light)
        next_state_C = press_C(current_number, current_light)
        
        if next_state_A and next_state_A not in visited:
            queue.append((next_state_A, sequence + ['A']))
        
        if next_state_B and next_state_B not in visited:
            queue.append((next_state_B, sequence + ['B']))
        
        if next_state_C and next_state_C not in visited:
            queue.append((next_state_C, sequence + ['C']))
    
    return None

# Find the shortest sequence
sequence = find_shortest_sequence(initial_state, target)
print(" → ".join(sequence))