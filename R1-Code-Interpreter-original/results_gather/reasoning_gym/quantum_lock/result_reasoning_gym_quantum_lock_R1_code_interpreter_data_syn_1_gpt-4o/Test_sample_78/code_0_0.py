from collections import deque

# Define the operations for each button
operations = {
    'A': lambda x: x - 3,
    'B': lambda x: x + 3,
    'C': lambda x: x + 2
}

# Initialize the queue with the starting state
queue = deque([(0, [])])  # (current number, sequence of buttons)

# Set to keep track of visited numbers to avoid cycles
visited = set()
visited.add(0)

# Target number
target = 19

# Perform BFS
while queue:
    current, sequence = queue.popleft()
    
    # Check if we have reached the target
    if current == target:
        print(" → ".join(sequence))
        break
    
    # Explore each button operation
    for button, operation in operations.items():
        next_number = operation(current)
        
        # If the next number is not visited, add it to the queue
        if next_number not in visited:
            visited.add(next_number)
            queue.append((next_number, sequence + [button]))