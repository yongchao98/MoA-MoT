from collections import deque

# Define the operations
def apply_operation(state, button):
    number, light = state
    if button == 'A':
        return (number - 3, not light)
    elif button == 'B':
        return (number + 2, not light)
    elif button == 'C':
        return (number * 3, not light)

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, [])])  # (state, sequence of operations)
    visited = set()
    
    while queue:
        current_state, sequence = queue.popleft()
        number, _ = current_state
        
        if number == target:
            return sequence
        
        for button in ['A', 'B', 'C']:
            new_state = apply_operation(current_state, button)
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, sequence + [button]))

# Initial state: (0, red)
start_state = (0, False)  # False for red, True for green
target_number = 39

# Find the shortest sequence
sequence = find_shortest_sequence(start_state, target_number)
print(" → ".join(sequence))