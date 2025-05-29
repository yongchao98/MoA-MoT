from collections import deque

# Define the initial state
start_number = 0
target_number = 30
initial_light = 'red'

# Define the operations
def press_A(number):
    return number + 3

def press_B(number):
    return number - 2

def press_C(number):
    return number * 2

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, initial_light, [])])  # (current number, light state, sequence of presses)
    visited = set()  # To keep track of visited states

    while queue:
        current_number, light, sequence = queue.popleft()

        # Check if we have reached the target
        if current_number == target:
            return sequence

        # If this state has been visited, skip it
        if (current_number, light) in visited:
            continue

        # Mark this state as visited
        visited.add((current_number, light))

        # Try pressing each button and add the new state to the queue
        # Press A
        new_number = press_A(current_number)
        new_light = 'green' if light == 'red' else 'red'
        queue.append((new_number, new_light, sequence + ['A']))

        # Press B
        new_number = press_B(current_number)
        new_light = 'green' if light == 'red' else 'red'
        queue.append((new_number, new_light, sequence + ['B']))

        # Press C (only if the light is green)
        if light == 'green':
            new_number = press_C(current_number)
            new_light = 'red'  # Light will toggle to red
            queue.append((new_number, new_light, sequence + ['C']))

# Find the shortest sequence
sequence = find_shortest_sequence(start_number, target_number)
print(' → '.join(sequence))