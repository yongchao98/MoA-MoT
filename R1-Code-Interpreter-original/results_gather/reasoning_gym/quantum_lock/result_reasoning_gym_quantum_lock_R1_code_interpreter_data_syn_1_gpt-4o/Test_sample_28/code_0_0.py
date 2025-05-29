from collections import deque

# Define the operations for each button
def button_A(n):
    return n - 2

def button_B(n):
    return n * 2

def button_C(n):
    return n + 3

# BFS to find the shortest sequence of button presses
def find_shortest_sequence(start, target):
    queue = deque([(start, "", 0)])  # (current_value, sequence_of_buttons, steps)
    visited = set()
    
    while queue:
        current_value, sequence, steps = queue.popleft()
        
        # If we reach the target, return the sequence
        if current_value == target:
            return sequence
        
        # If already visited this value, skip it
        if current_value in visited:
            continue
        
        # Mark this value as visited
        visited.add(current_value)
        
        # Explore the next states
        for button, operation in [('A', button_A), ('B', button_B), ('C', button_C)]:
            next_value = operation(current_value)
            next_sequence = sequence + (f" → {button}" if sequence else button)
            queue.append((next_value, next_sequence, steps + 1))
    
    return None

# Start the search
result = find_shortest_sequence(0, 16)
print(result)