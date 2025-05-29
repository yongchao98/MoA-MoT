from collections import deque

# Define the initial state
start_number = 0
target_number = 34
initial_light = 'red'  # Light starts as red

# Define the operations for each button
def press_A(number):
    return number + 3

def press_B(number):
    return number * 2

def press_C(number, light):
    if light == 'green':
        return number - 2
    return number

# Define the BFS function to find the shortest sequence
def find_shortest_sequence(start, target):
    # Queue for BFS: stores tuples of (current_number, light_state, sequence_of_buttons)
    queue = deque([(start, initial_light, [])])
    visited = set()  # To keep track of visited states

    while queue:
        current_number, light, sequence = queue.popleft()

        # Check if we have reached the target
        if current_number == target:
            return sequence

        # If this state has been visited, skip it
        if (current_number, light) in visited:
            continue
        visited.add((current_number, light))

        # Try pressing each button and add the new state to the queue
        # Press A
        new_number = press_A(current_number)
        new_light = 'green' if light == 'red' else 'red'
        queue.append((new_number, new_light, sequence + ['A']))

        # Press B
        new_number = press_B(current_number)
        new_light = 'green' if light == 'red' else 'red'
        queue.append((new_number, new_light, sequence + ['B']))

        # Press C (only if light is green)
        if light == 'green':
            new_number = press_C(current_number, light)
            new_light = 'red'  # Light toggles to red
            queue.append((new_number, new_light, sequence + ['C']))

# Find the shortest sequence
shortest_sequence = find_shortest_sequence(start_number, target_number)
print(" → ".join(shortest_sequence))