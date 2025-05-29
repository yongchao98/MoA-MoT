from collections import deque

def find_shortest_sequence(start, target):
    # Queue for BFS: stores tuples of (current number, light state, sequence of buttons)
    queue = deque([(start, 'red', [])])
    visited = set()  # To keep track of visited states

    while queue:
        current, light, sequence = queue.popleft()

        # If we reach the target, return the sequence
        if current == target:
            return sequence

        # If this state has been visited, skip it
        if (current, light) in visited:
            continue
        visited.add((current, light))

        # Calculate the next states for each button
        if light == 'green':
            # Button A: Add 2
            queue.append((current + 2, 'red', sequence + ['A']))
            # Button B: Subtract 3
            queue.append((current - 3, 'red', sequence + ['B']))
        
        # Button C: Multiply by 3 (works for both red and green)
        queue.append((current * 3, 'green' if light == 'red' else 'red', sequence + ['C']))

# Find the shortest sequence
sequence = find_shortest_sequence(0, 35)
print(" → ".join(sequence))