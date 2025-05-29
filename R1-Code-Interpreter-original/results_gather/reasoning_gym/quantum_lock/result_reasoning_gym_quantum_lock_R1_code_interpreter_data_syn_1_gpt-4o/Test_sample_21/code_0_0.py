from collections import deque

# Define the initial state
start = (0, 'red')  # (number, light state)
target = 30

# Define the button operations
def press_A(number, light):
    return number + 3, 'green' if light == 'red' else 'red'

def press_B(number, light):
    if light == 'green':
        return number - 3, 'red'
    return None

def press_C(number, light):
    if light == 'red':
        return number * 2, 'green'
    return None

# BFS to find the shortest sequence
queue = deque([(start, [])])  # (state, path)
visited = set()

while queue:
    (current_number, current_light), path = queue.popleft()
    
    if current_number == target:
        print(' → '.join(path))
        break
    
    # Try pressing each button
    for button, operation in [('A', press_A), ('B', press_B), ('C', press_C)]:
        result = operation(current_number, current_light)
        if result and result not in visited:
            visited.add(result)
            queue.append((result, path + [button]))