from collections import deque

# Define the initial state
start_number = 0
target_number = 20
initial_light = 'red'

# Define the operations for each button
def press_A(number, light):
    return number + 2, 'green' if light == 'red' else 'red'

def press_B(number, light):
    if light == 'green':
        return number - 3, 'red'
    return None

def press_C(number, light):
    if light == 'green':
        return number * 2, 'red'
    return None

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, initial_light, [])])  # (current number, current light, sequence of buttons)
    visited = set()

    while queue:
        current_number, current_light, sequence = queue.popleft()

        if current_number == target:
            return sequence

        if (current_number, current_light) in visited:
            continue

        visited.add((current_number, current_light))

        # Try pressing each button
        new_number, new_light = press_A(current_number, current_light)
        queue.append((new_number, new_light, sequence + ['A']))

        if current_light == 'green':
            result = press_B(current_number, current_light)
            if result:
                new_number, new_light = result
                queue.append((new_number, new_light, sequence + ['B']))

            result = press_C(current_number, current_light)
            if result:
                new_number, new_light = result
                queue.append((new_number, new_light, sequence + ['C']))

# Find the shortest sequence
sequence = find_shortest_sequence(start_number, target_number)
print(' → '.join(sequence))