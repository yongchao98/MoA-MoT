from collections import deque

# Define the button operations
operations = {
    'A': lambda x: x + 3,
    'B': lambda x: x + 2,
    'C': lambda x: x - 3
}

# Initialize the BFS queue
queue = deque([(0, [])])  # (current number, sequence of button presses)
visited = set()  # To keep track of visited numbers

# Perform BFS
while queue:
    current, sequence = queue.popleft()
    
    # Check if we've reached the target
    if current == 14:
        print(" → ".join(sequence))
        break
    
    # Explore the next states
    for button, operation in operations.items():
        next_number = operation(current)
        if next_number not in visited:
            visited.add(next_number)
            queue.append((next_number, sequence + [button]))