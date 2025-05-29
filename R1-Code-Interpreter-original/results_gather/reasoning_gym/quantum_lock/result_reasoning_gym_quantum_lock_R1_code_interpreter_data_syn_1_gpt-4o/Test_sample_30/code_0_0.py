from collections import deque

# Define the initial state
start = 0
target = 34
initial_light = 'red'

# Define the operations
def press_A(number, light):
    if light == 'red':
        return number + 2, 'green'
    return number, light

def press_B(number, light):
    return number + 3, 'red' if light == 'green' else 'green'

def press_C(number, light):
    return number * 3, 'red' if light == 'green' else 'green'

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, initial_light, [])])  # (current number, light state, sequence of operations)
    visited = set((start, initial_light))
    
    while queue:
        current, light, sequence = queue.popleft()
        
        if current == target:
            return sequence
        
        # Try each button
        for button, operation in [('A', press_A), ('B', press_B), ('C', press_C)]:
            next_number, next_light = operation(current, light)
            if (next_number, next_light) not in visited:
                visited.add((next_number, next_light))
                queue.append((next_number, next_light, sequence + [button]))
    
    return None

# Find and print the shortest sequence
sequence = find_shortest_sequence(start, target)
print(" → ".join(sequence))