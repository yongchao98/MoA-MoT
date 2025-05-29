from collections import deque

# Initial state
start_number = 0
target_number = 17
initial_light = 'red'  # Light starts as red

# Queue for BFS: (current_number, light_state, sequence_of_presses)
queue = deque([(start_number, initial_light, [])])

# Set to keep track of visited states
visited = set()
visited.add((start_number, initial_light))

# BFS to find the shortest sequence
while queue:
    current_number, light_state, sequence = queue.popleft()
    
    # Check if we've reached the target
    if current_number == target_number:
        print(" → ".join(sequence))
        break
    
    # Explore the next states
    # Button A: Subtract 2
    next_number = current_number - 2
    next_light = 'green' if light_state == 'red' else 'red'
    if (next_number, next_light) not in visited:
        visited.add((next_number, next_light))
        queue.append((next_number, next_light, sequence + ['A']))
    
    # Button B: Add 3 (only if light is green)
    if light_state == 'green':
        next_number = current_number + 3
        next_light = 'red'
        if (next_number, next_light) not in visited:
            visited.add((next_number, next_light))
            queue.append((next_number, next_light, sequence + ['B']))
    
    # Button C: Add 3
    next_number = current_number + 3
    next_light = 'green' if light_state == 'red' else 'red'
    if (next_number, next_light) not in visited:
        visited.add((next_number, next_light))
        queue.append((next_number, next_light, sequence + ['C']))