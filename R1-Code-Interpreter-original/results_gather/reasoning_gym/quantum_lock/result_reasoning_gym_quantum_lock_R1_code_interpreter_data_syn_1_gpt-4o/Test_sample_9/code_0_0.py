from collections import deque

# Define the operations
def apply_operation(state, button):
    if button == 'A':
        return state + 2
    elif button == 'B':
        return state - 3
    elif button == 'C':
        return state * 3

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, "", 0)])  # (current number, sequence of buttons, depth)
    visited = set()
    
    while queue:
        current, sequence, depth = queue.popleft()
        
        if current == target:
            return sequence
        
        if current in visited:
            continue
        
        visited.add(current)
        
        for button in ['A', 'B', 'C']:
            new_state = apply_operation(current, button)
            if new_state not in visited:
                queue.append((new_state, sequence + (button + " → "), depth + 1))
    
    return "No solution"

# Start the search
result = find_shortest_sequence(0, 19)
print(result.strip(" → "))  # Remove the trailing arrow