from collections import deque

# Define the initial state
initial_state = (0, 'red')  # (number, light color)
target_number = 15

# Define the operations for each button
def press_A(number):
    return number * 2

def press_B(number):
    return number + 1

def press_C(number):
    return number - 1

# BFS to find the shortest sequence of button presses
def find_shortest_sequence(initial_state, target_number):
    queue = deque([(initial_state, [])])  # (state, sequence of buttons)
    visited = set()
    
    while queue:
        (current_number, current_light), sequence = queue.popleft()
        
        # Check if we reached the target number
        if current_number == target_number:
            return sequence
        
        # If already visited this state, skip it
        if (current_number, current_light) in visited:
            continue
        visited.add((current_number, current_light))
        
        # Try pressing each button
        # Press A
        new_number = press_A(current_number)
        new_light = 'green' if current_light == 'red' else 'red'
        queue.append(((new_number, new_light), sequence + ['A']))
        
        # Press B
        new_number = press_B(current_number)
        new_light = 'green' if current_light == 'red' else 'red'
        queue.append(((new_number, new_light), sequence + ['B']))
        
        # Press C (only if light is green)
        if current_light == 'green':
            new_number = press_C(current_number)
            new_light = 'red'
            queue.append(((new_number, new_light), sequence + ['C']))

# Find the shortest sequence
sequence = find_shortest_sequence(initial_state, target_number)
print(" → ".join(sequence))