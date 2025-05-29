from collections import deque

# Define the initial state
initial_state = (0, 'red')  # (number, light)

# Define the target number
target = 27

# Define the button operations
def press_A(number):
    return number + 2

def press_B(number):
    return number * 3

def press_C(number):
    return number + 3

# Define the BFS queue
queue = deque([(initial_state, [])])  # (state, sequence of buttons)

# Define a set to keep track of visited states
visited = set()

# Perform BFS
while queue:
    (current_number, current_light), sequence = queue.popleft()
    
    # Check if we have reached the target
    if current_number == target:
        print('<<<' + ' → '.join(sequence) + '>>>')
        break
    
    # If the current state has been visited, skip it
    if (current_number, current_light) in visited:
        continue
    
    # Mark the current state as visited
    visited.add((current_number, current_light))
    
    # Generate new states based on the current light
    if current_light == 'red':
        # Press button A
        new_number = press_A(current_number)
        new_light = 'green'
        queue.append(((new_number, new_light), sequence + ['A']))
    else:
        # Press button B
        new_number = press_B(current_number)
        new_light = 'red'
        queue.append(((new_number, new_light), sequence + ['B']))
        
        # Press button C
        new_number = press_C(current_number)
        new_light = 'red'
        queue.append(((new_number, new_light), sequence + ['C']))