from collections import deque

def find_shortest_sequence(start, target):
    # Initial state
    initial_state = (start, 'red')
    # Queue for BFS: stores tuples of (current_number, light_color, sequence_of_buttons)
    queue = deque([(initial_state, [])])
    # Set to keep track of visited states
    visited = set()
    visited.add(initial_state)
    
    while queue:
        (current_number, light_color), sequence = queue.popleft()
        
        # Check if we've reached the target
        if current_number == target:
            return sequence
        
        # Generate possible next states
        # Button A: Add 1, toggle light
        next_number = current_number + 1
        next_light = 'green' if light_color == 'red' else 'red'
        next_state = (next_number, next_light)
        if next_state not in visited:
            visited.add(next_state)
            queue.append((next_state, sequence + ['A']))
        
        # Button B: Add 2, toggle light (only if light is green)
        if light_color == 'green':
            next_number = current_number + 2
            next_light = 'red'
            next_state = (next_number, next_light)
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, sequence + ['B']))
        
        # Button C: Multiply by 2, toggle light
        next_number = current_number * 2
        next_light = 'green' if light_color == 'red' else 'red'
        next_state = (next_number, next_light)
        if next_state not in visited:
            visited.add(next_state)
            queue.append((next_state, sequence + ['C']))

# Find the shortest sequence to reach the target
sequence = find_shortest_sequence(0, 16)
print(" → ".join(sequence))