from collections import deque

# Initial state
start_state = (0, 'red')
target_number = 9

# Define the button operations
def toggle(color):
    return 'green' if color == 'red' else 'red'

def press_A(state):
    number, color = state
    return (number - 1, toggle(color))

def press_B(state):
    number, color = state
    if color == 'red':
        return (number * 2, toggle(color))
    return None

def press_C(state):
    number, color = state
    return (number + 2, toggle(color))

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, [])])  # (state, path)
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        number, color = current_state
        
        if number == target:
            return path
        
        if current_state in visited:
            continue
        visited.add(current_state)
        
        # Try each button
        for button, operation in [('A', press_A), ('B', press_B), ('C', press_C)]:
            new_state = operation(current_state)
            if new_state and new_state not in visited:
                queue.append((new_state, path + [button]))
    
    return None

# Find the shortest sequence
sequence = find_shortest_sequence(start_state, target_number)
print(" → ".join(sequence))