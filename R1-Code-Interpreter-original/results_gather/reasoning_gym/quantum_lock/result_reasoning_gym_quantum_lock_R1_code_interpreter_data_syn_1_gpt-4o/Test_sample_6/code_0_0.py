from collections import deque

# Define the operations
def apply_operation(number, operation):
    if operation == 'A':
        return number + 3
    elif operation == 'B':
        return number * 2
    elif operation == 'C':
        return number - 2

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, [])])  # (current number, sequence of operations)
    visited = set()
    
    while queue:
        current, sequence = queue.popleft()
        
        if current == target:
            return sequence
        
        for operation in ['A', 'B', 'C']:
            next_number = apply_operation(current, operation)
            if next_number not in visited:
                visited.add(next_number)
                queue.append((next_number, sequence + [operation]))

# Start the search
start_number = 0
target_number = 20
sequence = find_shortest_sequence(start_number, target_number)

# Format the sequence for output
formatted_sequence = ' → '.join(sequence)
print(f"<<<{formatted_sequence}>>>")