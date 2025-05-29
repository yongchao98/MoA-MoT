from collections import deque

# Define the operations for each button
def press_A(n):
    return n - 2

def press_B(n):
    return n + 1

def press_C(n):
    return n * 2

# BFS to find the shortest sequence of button presses
def find_shortest_sequence(start, target):
    queue = deque([(start, "", 0)])  # (current_value, sequence_of_buttons, steps)
    visited = set()
    
    while queue:
        current_value, sequence, steps = queue.popleft()
        
        if current_value == target:
            return sequence
        
        if current_value in visited:
            continue
        
        visited.add(current_value)
        
        # Try each button press
        for button, operation in [('A', press_A), ('B', press_B), ('C', press_C)]:
            new_value = operation(current_value)
            new_sequence = sequence + (f" → {button}" if sequence else button)
            queue.append((new_value, new_sequence, steps + 1))
    
    return None

# Start the search
result = find_shortest_sequence(0, 11)
print(result)