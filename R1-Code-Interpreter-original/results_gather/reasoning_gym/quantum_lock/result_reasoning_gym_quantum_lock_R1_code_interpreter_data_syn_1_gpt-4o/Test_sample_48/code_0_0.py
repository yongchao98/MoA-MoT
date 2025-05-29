from collections import deque

# Define the initial state
start_number = 0
target_number = 14
initial_light = 'red'  # Light starts as red

# Define the button operations
def press_A(number, light):
    return number - 2, 'green' if light == 'red' else 'red'

def press_B(number, light):
    if light == 'green':
        return number - 3, 'red'
    return None

def press_C(number, light):
    return number + 2, 'green' if light == 'red' else 'red'

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, initial_light, [])])  # (current number, light color, sequence of buttons)
    visited = set()

    while queue:
        current_number, current_light, sequence = queue.popleft()

        # Check if we reached the target
        if current_number == target:
            return sequence

        # Mark the current state as visited
        visited.add((current_number, current_light))

        # Try pressing each button
        for button, operation in [('A', press_A), ('B', press_B), ('C', press_C)]:
            result = operation(current_number, current_light)
            if result:
                new_number, new_light = result
                if (new_number, new_light) not in visited:
                    queue.append((new_number, new_light, sequence + [button]))

# Find the shortest sequence
sequence = find_shortest_sequence(start_number, target_number)
print(" → ".join(sequence))