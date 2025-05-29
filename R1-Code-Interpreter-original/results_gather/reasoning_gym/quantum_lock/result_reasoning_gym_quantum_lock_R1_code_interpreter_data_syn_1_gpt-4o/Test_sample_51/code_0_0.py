from collections import deque

# Define the initial state
start_number = 0
target_number = 24
initial_light = 'red'

# Define the operations for each button
def press_A(number, light):
    return number + 3, 'green' if light == 'red' else 'red'

def press_B(number, light):
    return number - 2, 'green' if light == 'red' else 'red'

def press_C(number, light):
    if light == 'green':
        return number * 2, 'red'
    return number, light

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, initial_light, [])])  # (current number, light state, sequence of presses)
    visited = set()

    while queue:
        current_number, current_light, sequence = queue.popleft()

        # Check if we reached the target
        if current_number == target:
            return sequence

        # Avoid revisiting the same state
        if (current_number, current_light) in visited:
            continue
        visited.add((current_number, current_light))

        # Try each button
        for button, operation in [('A', press_A), ('B', press_B), ('C', press_C)]:
            new_number, new_light = operation(current_number, current_light)
            new_sequence = sequence + [button]
            queue.append((new_number, new_light, new_sequence))

# Find the shortest sequence
shortest_sequence = find_shortest_sequence(start_number, target_number)
print(" → ".join(shortest_sequence))