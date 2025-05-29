from collections import deque

def find_shortest_sequence():
    start = 0
    target = 8
    initial_light = 'red'
    
    # Queue will store tuples of (current_value, light_color, sequence_of_buttons)
    queue = deque([(start, initial_light, [])])
    visited = set((start, initial_light))
    
    while queue:
        current_value, light_color, sequence = queue.popleft()
        
        if current_value == target:
            return ' → '.join(sequence)
        
        # Try all possible button presses
        # Button A: Multiply by 2
        new_value = current_value * 2
        new_light = 'green' if light_color == 'red' else 'red'
        if (new_value, new_light) not in visited:
            visited.add((new_value, new_light))
            queue.append((new_value, new_light, sequence + ['A']))
        
        # Button B: Add 1 (only if light is green)
        if light_color == 'green':
            new_value = current_value + 1
            new_light = 'red'
            if (new_value, new_light) not in visited:
                visited.add((new_value, new_light))
                queue.append((new_value, new_light, sequence + ['B']))
        
        # Button C: Subtract 1
        new_value = current_value - 1
        new_light = 'green' if light_color == 'red' else 'red'
        if (new_value, new_light) not in visited:
            visited.add((new_value, new_light))
            queue.append((new_value, new_light, sequence + ['C']))
    
    return "No solution found"

print(find_shortest_sequence())