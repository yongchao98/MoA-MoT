from collections import deque

# Define the operations for each button
def button_A(n):
    return n * 2

def button_B(n):
    return n + 3

def button_C(n):
    return n + 2

# BFS to find the shortest sequence of button presses
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
        
        # Try each button and add the result to the queue
        queue.append((button_A(current), sequence + "A → ", depth + 1))
        queue.append((button_B(current), sequence + "B → ", depth + 1))
        queue.append((button_C(current), sequence + "C → ", depth + 1))
    
    return "No solution"

# Start the search
result = find_shortest_sequence(0, 27)
print(result.strip(" → "))  # Remove the trailing arrow