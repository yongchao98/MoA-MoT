from collections import deque

# Initial state
start_number = 0
start_light = 'red'
target_number = 8

# Define the operations for each button
def press_A(number, light):
    return number + 2, 'green' if light == 'red' else 'red'

def press_B(number, light):
    if light == 'green':
        return number * 2, 'red'
    return None

def press_C(number, light):
    if light == 'red':
        return number * 2, 'green'
    return None

# BFS to find the shortest sequence
def find_shortest_sequence(start_number, start_light, target_number):
    queue = deque([(start_number, start_light, [])])
    visited = set()

    while queue:
        number, light, path = queue.popleft()

        if number == target_number:
            return path

        if (number, light) in visited:
            continue

        visited.add((number, light))

        # Try pressing each button
        next_state = press_A(number, light)
        if next_state:
            queue.append((next_state[0], next_state[1], path + ['A']))

        next_state = press_B(number, light)
        if next_state:
            queue.append((next_state[0], next_state[1], path + ['B']))

        next_state = press_C(number, light)
        if next_state:
            queue.append((next_state[0], next_state[1], path + ['C']))

# Find the shortest sequence
sequence = find_shortest_sequence(start_number, start_light, target_number)
print(" → ".join(sequence))