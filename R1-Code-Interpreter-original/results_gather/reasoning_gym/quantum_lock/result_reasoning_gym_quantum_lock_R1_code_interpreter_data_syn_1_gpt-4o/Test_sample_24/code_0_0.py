from collections import deque

# Define the operations for each button
def apply_operation(state, button):
    number, _ = state
    if button == 'A':
        return number + 2
    elif button == 'B':
        return number * 2
    elif button == 'C':
        return number * 3

# BFS to find the shortest sequence of button presses
def find_shortest_sequence(start, target):
    queue = deque([(start, [])])  # (current state, path)
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        current_number, _ = current_state
        
        if current_number == target:
            return path
        
        for button in ['A', 'B', 'C']:
            new_number = apply_operation(current_state, button)
            new_state = (new_number, not current_state[1])  # Toggle light state
            
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + [button]))

# Initial state: (0, red)
start_state = (0, False)  # False for red, True for green
target_number = 28

# Find the shortest sequence
sequence = find_shortest_sequence(start_state, target_number)
print(" → ".join(sequence))