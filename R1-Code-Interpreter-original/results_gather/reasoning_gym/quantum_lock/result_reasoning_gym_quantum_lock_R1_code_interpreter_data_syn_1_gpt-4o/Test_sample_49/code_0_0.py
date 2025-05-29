from collections import deque

# Define the initial state
start_number = 0
target_number = 26
initial_light = 'red'  # 'red' or 'green'

# Define the button operations
def press_button(state, button):
    number, light = state
    if button == 'A':
        return (number - 3, 'green' if light == 'red' else 'red')
    elif button == 'B' and light == 'green':
        return (number - 3, 'red')
    elif button == 'C':
        return (number + 2, 'green' if light == 'red' else 'red')
    return None

# BFS to find the shortest sequence
def find_shortest_sequence(start, target):
    queue = deque([(start, [])])  # (state, sequence of buttons)
    visited = set()
    
    while queue:
        current_state, sequence = queue.popleft()
        number, light = current_state
        
        if number == target:
            return sequence
        
        if current_state in visited:
            continue
        visited.add(current_state)
        
        for button in ['A', 'B', 'C']:
            new_state = press_button(current_state, button)
            if new_state and new_state not in visited:
                queue.append((new_state, sequence + [button]))
    
    return None

# Initial state
initial_state = (start_number, initial_light)

# Find the shortest sequence
shortest_sequence = find_shortest_sequence(initial_state, target_number)

# Format the output
if shortest_sequence:
    result = ' → '.join(shortest_sequence)
else:
    result = "No solution found"

print(f"<<<{result}>>>")