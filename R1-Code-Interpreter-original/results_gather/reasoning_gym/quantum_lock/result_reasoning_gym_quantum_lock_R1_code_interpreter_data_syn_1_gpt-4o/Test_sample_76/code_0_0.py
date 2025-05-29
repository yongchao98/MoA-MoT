from collections import deque

# Define the initial state
start = (0, 'red')  # (number, light color)
target = 13

# Define the button operations
def press_A(number, light):
    return (number * 2, 'green' if light == 'red' else 'red')

def press_B(number, light):
    if light == 'red':
        return (number + 1, 'green')
    return None

def press_C(number, light):
    return (number - 1, 'green' if light == 'red' else 'red')

# BFS to find the shortest sequence
queue = deque([(start, [])])  # (state, path)
visited = set()

while queue:
    (current_number, current_light), path = queue.popleft()
    
    if current_number == target:
        print(' → '.join(path))
        break
    
    # Try pressing each button
    # Press A
    new_state = press_A(current_number, current_light)
    if new_state not in visited:
        visited.add(new_state)
        queue.append((new_state, path + ['A']))
    
    # Press B
    if current_light == 'red':
        new_state = press_B(current_number, current_light)
        if new_state not in visited:
            visited.add(new_state)
            queue.append((new_state, path + ['B']))
    
    # Press C
    new_state = press_C(current_number, current_light)
    if new_state not in visited:
        visited.add(new_state)
        queue.append((new_state, path + ['C']))