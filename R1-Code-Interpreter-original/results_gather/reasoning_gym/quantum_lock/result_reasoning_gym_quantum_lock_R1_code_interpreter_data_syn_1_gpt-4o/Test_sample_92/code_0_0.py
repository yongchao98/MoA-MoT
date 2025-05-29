from collections import deque

# Define the operations for each button
def operation_A(n):
    return n * 2

def operation_B(n):
    return n + 3

def operation_C(n):
    return n - 2

# BFS to find the shortest sequence of operations
def find_shortest_sequence(start, target):
    # Queue for BFS: stores tuples of (current number, sequence of operations)
    queue = deque([(start, [])])
    visited = set()  # To keep track of visited numbers

    while queue:
        current, sequence = queue.popleft()

        # If we reach the target, return the sequence
        if current == target:
            return sequence

        # If the current number has been visited, skip it
        if current in visited:
            continue

        # Mark the current number as visited
        visited.add(current)

        # Explore the next states
        next_states = [
            (operation_A(current), sequence + ['A']),
            (operation_B(current), sequence + ['B']),
            (operation_C(current), sequence + ['C'])
        ]

        for next_state, next_sequence in next_states:
            if next_state not in visited:
                queue.append((next_state, next_sequence))

# Find the shortest sequence from 0 to 19
sequence = find_shortest_sequence(0, 19)
print(" → ".join(sequence))