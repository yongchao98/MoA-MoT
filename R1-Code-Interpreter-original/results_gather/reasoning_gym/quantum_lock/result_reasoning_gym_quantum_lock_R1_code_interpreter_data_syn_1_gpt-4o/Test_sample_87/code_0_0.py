from collections import deque

# Define the initial state
start_number = 0
target_number = 13
initial_light = 'red'

# Define the button operations
def press_A(number, light):
    return number + 3, 'green' if light == 'red' else 'red'

def press_B(number, light):
    if light == 'red':
        return number - 2, 'green'
    return number, light

def press_C(number, light):
    return number - 3, 'green' if light == 'red' else 'red'

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, initial_light, [])])  # (current number, current light, sequence of presses)
    visited = set()  # To keep track of visited states

    while queue:
        current_number, current_light, sequence = queue.popleft()

        # Check if we reached the target
        if current_number == target:
            return sequence

        # Check if this state has been visited
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
shortest_sequence = find_shortest_sequence(start_number, target_number)
print(' → '.join(shortest_sequence))