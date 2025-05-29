from collections import deque

# Initial state
start = 0
target = 36
initial_light = 'red'

# Define the operations for each button
def press_A(number, light):
    if light == 'green':
        return number - 3, 'red'
    return number, light

def press_B(number, light):
    if light == 'green':
        return number * 3, 'red'
    return number, light

def press_C(number, light):
    return number + 3, 'green' if light == 'red' else 'red'

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, initial_light, [])])
    visited = set()

    while queue:
        current_number, current_light, sequence = queue.popleft()

        if current_number == target:
            return sequence

        if (current_number, current_light) in visited:
            continue

        visited.add((current_number, current_light))

        # Try each button
        new_number, new_light = press_A(current_number, current_light)
        queue.append((new_number, new_light, sequence + ['A']))

        new_number, new_light = press_B(current_number, current_light)
        queue.append((new_number, new_light, sequence + ['B']))

        new_number, new_light = press_C(current_number, current_light)
        queue.append((new_number, new_light, sequence + ['C']))

# Find the shortest sequence
sequence = find_shortest_sequence(start, target)
print(" → ".join(sequence))