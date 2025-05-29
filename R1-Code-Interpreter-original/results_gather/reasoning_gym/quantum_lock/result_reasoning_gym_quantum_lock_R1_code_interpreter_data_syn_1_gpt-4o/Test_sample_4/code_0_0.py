from collections import deque

# Initial state
start_number = 0
start_light = 'red'
target_number = 22

# Define the operations for each button
def apply_button(number, light, button):
    if button == 'A' and light == 'red':
        return number + 2, 'green'
    elif button == 'B':
        return number - 2, 'green' if light == 'red' else 'red'
    elif button == 'C':
        return number * 2, 'green' if light == 'red' else 'red'
    return None

# BFS to find the shortest sequence
def find_shortest_sequence(start_number, start_light, target_number):
    queue = deque([(start_number, start_light, [])])
    visited = set()

    while queue:
        number, light, sequence = queue.popleft()

        # Check if we reached the target
        if number == target_number:
            return sequence

        # Explore each button
        for button in ['A', 'B', 'C']:
            new_state = apply_button(number, light, button)
            if new_state:
                new_number, new_light = new_state
                if (new_number, new_light) not in visited:
                    visited.add((new_number, new_light))
                    queue.append((new_number, new_light, sequence + [button]))

# Find the sequence
sequence = find_shortest_sequence(start_number, start_light, target_number)
print(" → ".join(sequence))