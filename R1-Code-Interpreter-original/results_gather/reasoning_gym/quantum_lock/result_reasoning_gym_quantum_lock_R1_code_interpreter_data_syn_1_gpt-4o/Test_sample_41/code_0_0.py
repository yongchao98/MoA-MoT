from collections import deque

# Define the operations for each button
def button_A(n):
    return n * 2

def button_B(n):
    return n + 3

def button_C(n):
    return n - 2

# BFS to find the shortest sequence of button presses
def find_shortest_sequence(start, target):
    # Queue for BFS: stores tuples of (current number, sequence of buttons)
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

        # Explore the next states by pressing each button
        next_states = [
            (button_A(current), sequence + ['A']),
            (button_B(current), sequence + ['B']),
            (button_C(current), sequence + ['C'])
        ]

        # Add valid next states to the queue
        for next_state, next_sequence in next_states:
            if next_state not in visited:
                queue.append((next_state, next_sequence))

# Find the shortest sequence from 0 to 17
sequence = find_shortest_sequence(0, 17)
print(" → ".join(sequence))