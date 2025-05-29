from collections import deque

# Initial state
start_number = 0
start_light = 'red'
target_number = 26

# Queue for BFS: (current_number, current_light, sequence_of_buttons)
queue = deque([(start_number, start_light, [])])

# Visited states to avoid revisiting
visited = set()
visited.add((start_number, start_light))

# BFS loop
while queue:
    current_number, current_light, sequence = queue.popleft()
    
    # Check if we reached the target
    if current_number == target_number:
        print(" → ".join(sequence))
        break
    
    # Apply button A (Multiply by 3)
    new_number = current_number * 3
    new_light = 'green' if current_light == 'red' else 'red'
    if (new_number, new_light) not in visited:
        visited.add((new_number, new_light))
        queue.append((new_number, new_light, sequence + ['A']))
    
    # Apply button B (Add 2 if red)
    if current_light == 'red':
        new_number = current_number + 2
        new_light = 'green'
        if (new_number, new_light) not in visited:
            visited.add((new_number, new_light))
            queue.append((new_number, new_light, sequence + ['B']))
    
    # Apply button C (Add 2 if green)
    if current_light == 'green':
        new_number = current_number + 2
        new_light = 'red'
        if (new_number, new_light) not in visited:
            visited.add((new_number, new_light))
            queue.append((new_number, new_light, sequence + ['C']))