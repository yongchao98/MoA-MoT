from collections import deque

# Define the operations
def apply_operation(value, operation):
    if operation == 'A':
        return value + 2
    elif operation == 'B':
        return value - 1
    elif operation == 'C':
        return value * 2

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, "", "red")])  # (current_value, sequence, light_state)
    visited = set()
    
    while queue:
        current_value, sequence, light_state = queue.popleft()
        
        if current_value == target:
            return sequence
        
        if current_value < 0 or current_value > target * 2:
            continue
        
        if (current_value, light_state) in visited:
            continue
        
        visited.add((current_value, light_state))
        
        for operation in ['A', 'B', 'C']:
            new_value = apply_operation(current_value, operation)
            new_sequence = sequence + (operation + " → " if sequence else operation)
            new_light_state = "green" if light_state == "red" else "red"
            queue.append((new_value, new_sequence, new_light_state))
    
    return None

# Find the shortest sequence
result = find_shortest_sequence(0, 16)
print(result)