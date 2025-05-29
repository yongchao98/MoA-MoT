from collections import deque

def solve_jug_problem():
    # Initial state: all jugs are empty
    initial_state = (0, 0, 0)
    # Target amount of water
    target = 3
    # Jug capacities
    capacities = (12, 6, 13)
    
    # Queue for BFS: stores (state, path)
    queue = deque([(initial_state, [])])
    # Set to keep track of visited states
    visited = set()
    visited.add(initial_state)
    
    while queue:
        (a, b, c), path = queue.popleft()
        
        # Check if we have reached the target
        if a == target or b == target or c == target:
            return path
        
        # Generate all possible next states
        next_states = []
        
        # Fill operations
        next_states.append(((capacities[0], b, c), path + ['fill A']))
        next_states.append(((a, capacities[1], c), path + ['fill B']))
        next_states.append(((a, b, capacities[2]), path + ['fill C']))
        
        # Empty operations
        next_states.append(((0, b, c), path + ['empty A']))
        next_states.append(((a, 0, c), path + ['empty B']))
        next_states.append(((a, b, 0), path + ['empty C']))
        
        # Pour operations
        # Pour A -> B
        pour = min(a, capacities[1] - b)
        next_states.append(((a - pour, b + pour, c), path + ['pour A->B']))
        
        # Pour A -> C
        pour = min(a, capacities[2] - c)
        next_states.append(((a - pour, b, c + pour), path + ['pour A->C']))
        
        # Pour B -> A
        pour = min(b, capacities[0] - a)
        next_states.append(((a + pour, b - pour, c), path + ['pour B->A']))
        
        # Pour B -> C
        pour = min(b, capacities[2] - c)
        next_states.append(((a, b - pour, c + pour), path + ['pour B->C']))
        
        # Pour C -> A
        pour = min(c, capacities[0] - a)
        next_states.append(((a + pour, b, c - pour), path + ['pour C->A']))
        
        # Pour C -> B
        pour = min(c, capacities[1] - b)
        next_states.append(((a, b + pour, c - pour), path + ['pour C->B']))
        
        # Add valid next states to the queue
        for state, new_path in next_states:
            if state not in visited:
                visited.add(state)
                queue.append((state, new_path))

# Execute the function and print the result
result = solve_jug_problem()
print(result)