from collections import deque

# Define the initial state
start = 0
target = 8

# Define the operations for each button
operations = {
    'A': lambda x: x * 2,
    'B': lambda x: x + 1,
    'C': lambda x: x - 2
}

# Use a queue to perform BFS
queue = deque([(start, [])])  # (current number, sequence of buttons)
visited = set()

while queue:
    current, sequence = queue.popleft()
    
    # If we reach the target, return the sequence
    if current == target:
        print(" → ".join(sequence))
        break
    
    # Explore each button operation
    for button, operation in operations.items():
        next_number = operation(current)
        
        # Check if the next number is valid and not visited
        if next_number not in visited and next_number >= 0:
            visited.add(next_number)
            queue.append((next_number, sequence + [button]))