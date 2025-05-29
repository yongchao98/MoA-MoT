from collections import deque

# Initial state
start_number = 0
start_light = 'red'
target_number = 32

# Define the operations for each button
def press_A(number):
    return number * 2

def press_B(number):
    return number + 2

def press_C(number):
    return number * 2

# BFS to find the shortest sequence
def find_shortest_sequence(start_number, start_light, target_number):
    # Queue for BFS: (current_number, current_light, sequence_of_buttons)
    queue = deque([(start_number, start_light, [])])
    visited = set()  # To avoid revisiting the same state

    while queue:
        current_number, current_light, sequence = queue.popleft()

        # Check if we reached the target
        if current_number == target_number:
            return sequence

        # If the current state is already visited, skip it
        if (current_number, current_light) in visited:
            continue
        visited.add((current_number, current_light))

        # Explore the next states
        if current_light == 'red':
            # Press A
            next_number = press_A(current_number)
            next_light = 'green'
            queue.append((next_number, next_light, sequence + ['A']))
        else:
            # Press B
            next_number = press_B(current_number)
            next_light = 'red'
            queue.append((next_number, next_light, sequence + ['B']))

            # Press C
            next_number = press_C(current_number)
            next_light = 'red'
            queue.append((next_number, next_light, sequence + ['C']))

# Find the shortest sequence
shortest_sequence = find_shortest_sequence(start_number, start_light, target_number)
print(" → ".join(shortest_sequence))