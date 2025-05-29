from collections import deque

# Define the initial state
initial_state = (0, 'red')  # (number, light color)
target_number = 26

# Define the button operations
def press_A(number, light):
    return number + 2, 'green' if light == 'red' else 'red'

def press_B(number, light):
    if light == 'red':
        return number * 2, 'green'
    return None

def press_C(number, light):
    if light == 'green':
        return number * 2, 'red'
    return None

# BFS to find the shortest sequence
def find_shortest_sequence():
    queue = deque([(initial_state, [])])  # (state, sequence of buttons)
    visited = set()
    visited.add(initial_state)

    while queue:
        (current_number, current_light), sequence = queue.popleft()

        # Check if we reached the target
        if current_number == target_number:
            return sequence

        # Try pressing each button
        for button, operation in [('A', press_A), ('B', press_B), ('C', press_C)]:
            new_state = operation(current_number, current_light)
            if new_state and new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, sequence + [button]))

# Find the shortest sequence
shortest_sequence = find_shortest_sequence()
print(" → ".join(shortest_sequence))