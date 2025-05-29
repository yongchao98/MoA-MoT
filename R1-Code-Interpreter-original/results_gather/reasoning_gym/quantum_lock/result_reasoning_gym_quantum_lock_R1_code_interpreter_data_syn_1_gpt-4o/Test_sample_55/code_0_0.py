from collections import deque

def find_shortest_sequence(start, target):
    # Initial state: (current number, light color, sequence of presses)
    initial_state = (start, 'red', [])
    queue = deque([initial_state])
    visited = set()

    while queue:
        current_number, light_color, sequence = queue.popleft()

        # Check if we reached the target
        if current_number == target:
            return sequence

        # Generate possible next states
        # Button A: Subtract 1
        next_number = current_number - 1
        next_light_color = 'green' if light_color == 'red' else 'red'
        if (next_number, next_light_color) not in visited:
            visited.add((next_number, next_light_color))
            queue.append((next_number, next_light_color, sequence + ['A']))

        # Button B: Add 2
        next_number = current_number + 2
        next_light_color = 'green' if light_color == 'red' else 'red'
        if (next_number, next_light_color) not in visited:
            visited.add((next_number, next_light_color))
            queue.append((next_number, next_light_color, sequence + ['B']))

        # Button C: Multiply by 2 (only if light is green)
        if light_color == 'green':
            next_number = current_number * 2
            next_light_color = 'red'  # Light toggles to red
            if (next_number, next_light_color) not in visited:
                visited.add((next_number, next_light_color))
                queue.append((next_number, next_light_color, sequence + ['C']))

# Find the shortest sequence to reach the target
sequence = find_shortest_sequence(0, 9)
print(" → ".join(sequence))