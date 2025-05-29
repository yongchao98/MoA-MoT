from collections import deque

# Define the operations
def apply_operation(value, button):
    if button == 'A':
        return value - 2
    elif button == 'B':
        return value + 3
    elif button == 'C':
        return value * 2

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, [])])  # (current_value, sequence_of_buttons)
    visited = set()
    
    while queue:
        current_value, sequence = queue.popleft()
        
        if current_value == target:
            return sequence
        
        if current_value in visited:
            continue
        
        visited.add(current_value)
        
        for button in ['A', 'B', 'C']:
            new_value = apply_operation(current_value, button)
            if new_value not in visited:
                queue.append((new_value, sequence + [button]))
    
    return None

# Find the shortest sequence from 0 to 13
sequence = find_shortest_sequence(0, 13)
print(" → ".join(sequence))