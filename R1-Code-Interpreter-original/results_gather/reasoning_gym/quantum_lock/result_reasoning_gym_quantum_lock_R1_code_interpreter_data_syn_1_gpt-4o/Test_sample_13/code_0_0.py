from collections import deque

# Define the initial state
start = 0
target = 18
initial_light = 'red'

# Define the operations for each button
def press_A(number):
    return number * 2

def press_B(number):
    return number - 2

def press_C(number):
    return number + 2

# BFS to find the shortest sequence of button presses
def find_shortest_sequence(start, target):
    # Queue for BFS: stores tuples of (current_number, current_light, sequence_of_buttons)
    queue = deque([(start, initial_light, [])])
    visited = set()  # To keep track of visited states

    while queue:
        current_number, current_light, sequence = queue.popleft()

        # Check if we have reached the target
        if current_number == target:
            return sequence

        # If this state has been visited, skip it
        if (current_number, current_light) in visited:
            continue

        # Mark this state as visited
        visited.add((current_number, current_light))

        # Try pressing each button and enqueue the resulting state
        # Press A
        new_number = press_A(current_number)
        new_light = 'green' if current_light == 'red' else 'red'
        queue.append((new_number, new_light, sequence + ['A']))

        # Press B (only if the light is green)
        if current_light == 'green':
            new_number = press_B(current_number)
            new_light = 'red'  # Light toggles
            queue.append((new_number, new_light, sequence + ['B']))

        # Press C
        new_number = press_C(current_number)
        new_light = 'green' if current_light == 'red' else 'red'
        queue.append((new_number, new_light, sequence + ['C']))

# Find the shortest sequence
shortest_sequence = find_shortest_sequence(start, target)
print(" → ".join(shortest_sequence))