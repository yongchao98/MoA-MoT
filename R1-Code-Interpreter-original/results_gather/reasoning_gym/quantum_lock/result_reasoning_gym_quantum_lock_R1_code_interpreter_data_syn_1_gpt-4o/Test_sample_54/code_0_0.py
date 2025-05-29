from collections import deque

# Define the operations
def operation_A(n):
    return n + 2

def operation_B(n):
    return n * 2

def operation_C(n):
    return n + 1

# BFS to find the shortest path
def find_shortest_sequence(start, target):
    queue = deque([(start, "", 0)])  # (current_value, sequence, steps)
    visited = set()
    
    while queue:
        current_value, sequence, steps = queue.popleft()
        
        if current_value == target:
            return sequence
        
        if current_value not in visited:
            visited.add(current_value)
            
            # Apply operations
            queue.append((operation_A(current_value), sequence + "A → ", steps + 1))
            queue.append((operation_B(current_value), sequence + "B → ", steps + 1))
            queue.append((operation_C(current_value), sequence + "C → ", steps + 1))

# Start from 0 and find the sequence to reach 8
result = find_shortest_sequence(0, 8)
print(result.strip(" → "))