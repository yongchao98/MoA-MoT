from collections import deque

# Define the initial state
start = 0
target = 8
initial_light = 'red'

# Define the operations for each button
def press_A(number):
    return number + 1

def press_B(number, light):
    if light == 'red':
        return number - 2
    return None

def press_C(number):
    return number * 2

# Define the BFS function to find the shortest sequence
def find_shortest_sequence(start, target, initial_light):
    # Queue for BFS: stores tuples of (current_number, current_light, sequence_of_buttons)
    queue = deque([(start, initial_light, [])])
    visited = set()  # To keep track of visited states

    while queue:
        current_number, current_light, sequence = queue.popleft()

        # Check if we have reached the target
        if current_number == target:
            return sequence

        # Toggle the light
        next_light = 'green' if current_light == 'red' else 'red'

        # Try pressing each button and add the new state to the queue
        # Press A
        new_number = press_A(current_number)
        if (new_number, next_light) not in visited:
            visited.add((new_number, next_light))
            queue.append((new_number, next_light, sequence + ['A']))

        # Press B
        new_number = press_B(current_number, current_light)
        if new_number is not None and (new_number, next_light) not in visited:
            visited.add((new_number, next_light))
            queue.append((new_number, next_light, sequence + ['B']))

        # Press C
        new_number = press_C(current_number)
        if (new_number, next_light) not in visited:
            visited.add((new_number, next_light))
            queue.append((new_number, next_light, sequence + ['C']))

# Find the shortest sequence
sequence = find_shortest_sequence(start, target, initial_light)
print(" → ".join(sequence))