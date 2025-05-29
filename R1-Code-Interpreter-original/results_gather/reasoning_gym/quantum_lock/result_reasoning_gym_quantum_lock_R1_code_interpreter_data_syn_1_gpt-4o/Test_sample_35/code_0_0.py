from collections import deque

# Define the initial state
start_number = 0
start_light = 'red'
target_number = 8

# Define the button operations
def press_A(number, light):
    return number * 2, 'green' if light == 'red' else 'red'

def press_B(number, light):
    if light == 'green':
        return number + 1, 'green' if light == 'red' else 'red'
    return None

def press_C(number, light):
    return number - 1, 'green' if light == 'red' else 'red'

# BFS to find the shortest sequence
def find_shortest_sequence(start_number, start_light, target_number):
    queue = deque([(start_number, start_light, [])])
    visited = set()

    while queue:
        current_number, current_light, path = queue.popleft()

        if current_number == target_number:
            return path

        if (current_number, current_light) in visited:
            continue

        visited.add((current_number, current_light))

        # Try pressing each button
        new_number, new_light = press_A(current_number, current_light)
        queue.append((new_number, new_light, path + ['A']))

        if current_light == 'green':
            result = press_B(current_number, current_light)
            if result:
                new_number, new_light = result
                queue.append((new_number, new_light, path + ['B']))

        new_number, new_light = press_C(current_number, current_light)
        queue.append((new_number, new_light, path + ['C']))

# Find and print the shortest sequence
shortest_sequence = find_shortest_sequence(start_number, start_light, target_number)
print(' → '.join(shortest_sequence))