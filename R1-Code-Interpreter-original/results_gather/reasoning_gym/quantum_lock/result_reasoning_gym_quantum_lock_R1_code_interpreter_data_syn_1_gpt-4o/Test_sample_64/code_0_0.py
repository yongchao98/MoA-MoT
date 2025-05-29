from collections import deque

# Define the initial state
initial_state = (0, 'red')  # (number, light color)
target = 21

# Define the button operations
def press_A(number, light):
    if light == 'red':
        return number - 3, 'green'
    return number, light

def press_B(number, light):
    if light == 'green':
        return number + 3, 'red'
    return number, light

def press_C(number, light):
    return number * 2, 'green' if light == 'red' else 'red'

# BFS to find the shortest sequence
def find_shortest_sequence(initial_state, target):
    queue = deque([(initial_state, [])])  # (state, path)
    visited = set()

    while queue:
        (current_number, current_light), path = queue.popleft()

        # Check if we reached the target
        if current_number == target:
            return path

        # If already visited this state, skip it
        if (current_number, current_light) in visited:
            continue
        visited.add((current_number, current_light))

        # Explore the next states
        next_states = [
            (press_A(current_number, current_light), path + ['A']),
            (press_B(current_number, current_light), path + ['B']),
            (press_C(current_number, current_light), path + ['C'])
        ]

        for state, new_path in next_states:
            queue.append((state, new_path))

# Find the shortest sequence
shortest_sequence = find_shortest_sequence(initial_state, target)
print(" → ".join(shortest_sequence))