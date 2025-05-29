def find_sequence(start, target):
    from collections import deque

    # Each state is represented as (current_value, light_state, sequence_of_buttons)
    queue = deque([(start, 'red', [])])
    visited = set()

    while queue:
        current_value, light_state, sequence = queue.popleft()

        # If we reach the target, return the sequence
        if current_value == target:
            return sequence

        # If this state has been visited, skip it
        if (current_value, light_state) in visited:
            continue
        visited.add((current_value, light_state))

        # Generate possible next states
        if light_state == 'red':
            # Button A: Add 2, toggle to green
            queue.append((current_value + 2, 'green', sequence + ['A']))
        else:
            # Button C: Add 2, toggle to red
            queue.append((current_value + 2, 'red', sequence + ['C']))

        # Button B: Multiply by 2, toggle light
        new_light_state = 'green' if light_state == 'red' else 'red'
        queue.append((current_value * 2, new_light_state, sequence + ['B']))

# Find the sequence
sequence = find_sequence(0, 30)
print(" → ".join(sequence))