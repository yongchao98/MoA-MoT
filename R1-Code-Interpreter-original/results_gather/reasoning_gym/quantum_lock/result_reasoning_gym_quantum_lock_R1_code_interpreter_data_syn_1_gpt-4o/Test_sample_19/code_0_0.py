from collections import deque

def find_shortest_sequence(start, target):
    # Initial state
    initial_state = (start, 'red')
    # Queue for BFS
    queue = deque([(initial_state, [])])
    # Visited set to avoid revisiting states
    visited = set()
    visited.add(initial_state)
    
    while queue:
        (current_number, light_color), path = queue.popleft()
        
        # Check if we reached the target
        if current_number == target:
            return path
        
        # Generate possible next states
        # Button A: Subtract 1
        next_number = current_number - 1
        next_light = 'green' if light_color == 'red' else 'red'
        next_state = (next_number, next_light)
        if next_state not in visited:
            visited.add(next_state)
            queue.append((next_state, path + ['A']))
        
        # Button B: Add 1 (only if light is green)
        if light_color == 'green':
            next_number = current_number + 1
            next_light = 'red'
            next_state = (next_number, next_light)
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + ['B']))
        
        # Button C: Multiply by 2
        next_number = current_number * 2
        next_light = 'green' if light_color == 'red' else 'red'
        next_state = (next_number, next_light)
        if next_state not in visited:
            visited.add(next_state)
            queue.append((next_state, path + ['C']))

# Find the shortest sequence
sequence = find_shortest_sequence(0, 14)
print(" → ".join(sequence))