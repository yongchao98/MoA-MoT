from collections import deque

def find_shortest_sequence(start, target):
    # Queue for BFS: stores tuples of (current number, light state, sequence of buttons)
    queue = deque([(start, 'red', [])])
    visited = set()  # To avoid revisiting the same state

    while queue:
        current, light, sequence = queue.popleft()

        # If we reach the target, return the sequence
        if current == target:
            return sequence

        # If this state has been visited, skip it
        if (current, light) in visited:
            continue
        visited.add((current, light))

        # Determine the next light state
        next_light = 'green' if light == 'red' else 'red'

        # Try each button and add the resulting state to the queue
        # Button A: Subtract 2 (when green)
        if light == 'green':
            queue.append((current - 2, next_light, sequence + ['A']))

        # Button B: Add 3 (when any)
        queue.append((current + 3, next_light, sequence + ['B']))

        # Button C: Subtract 3 (when green)
        if light == 'green':
            queue.append((current - 3, next_light, sequence + ['C']))

# Find the shortest sequence
sequence = find_shortest_sequence(0, 29)
print(" → ".join(sequence))