def find_shortest_sequence(start, target):
    from collections import deque

    # Each state is represented as (current_value, sequence_of_buttons)
    queue = deque([(start, [])])
    visited = set()

    while queue:
        current_value, sequence = queue.popleft()

        # If we reach the target, return the sequence
        if current_value == target:
            return sequence

        # If already visited this state, skip it
        if current_value in visited:
            continue

        # Mark this state as visited
        visited.add(current_value)

        # Try each button and add the new state to the queue
        for button, operation in [('A', -3), ('B', -2), ('C', 3)]:
            new_value = current_value + operation
            if new_value not in visited:
                queue.append((new_value, sequence + [button]))

# Find the shortest sequence from 0 to 23
sequence = find_shortest_sequence(0, 23)
print(" → ".join(sequence))