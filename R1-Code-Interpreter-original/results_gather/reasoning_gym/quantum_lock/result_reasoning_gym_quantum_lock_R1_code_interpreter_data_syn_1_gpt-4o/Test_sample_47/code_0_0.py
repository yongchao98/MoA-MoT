from collections import deque

# Define the initial state
start = 0
target = 10
initial_light = 'red'

# Define the operations for each button
def press_A(number, light):
    if light == 'red':
        return number * 2, 'green'
    return number, light

def press_B(number, light):
    if light == 'green':
        return number - 1, 'red'
    return number, light

def press_C(number, light):
    return number + 1, 'green' if light == 'red' else 'red'

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, initial_light, [])])  # (current number, light state, sequence of buttons)
    visited = set()

    while queue:
        current, light, sequence = queue.popleft()

        # Check if we reached the target
        if current == target:
            return sequence

        # Avoid revisiting the same state
        if (current, light) in visited:
            continue
        visited.add((current, light))

        # Try each button
        next_number, next_light = press_A(current, light)
        queue.append((next_number, next_light, sequence + ['A']))

        next_number, next_light = press_B(current, light)
        queue.append((next_number, next_light, sequence + ['B']))

        next_number, next_light = press_C(current, light)
        queue.append((next_number, next_light, sequence + ['C']))

# Find the shortest sequence
sequence = find_shortest_sequence(start, target)
print(" → ".join(sequence))