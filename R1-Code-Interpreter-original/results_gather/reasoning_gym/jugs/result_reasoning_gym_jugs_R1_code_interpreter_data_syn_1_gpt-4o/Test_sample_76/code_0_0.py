from collections import deque

def solve_jug_problem():
    # Initial state: all jugs are empty
    initial_state = (0, 0, 0)
    # Target amount of water
    target = 4
    # Jug capacities
    capacities = (8, 8, 13)
    
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
        if a < capacities[0]:
            next_states.append(((capacities[0], b, c), path + ['fill A']))
        if b < capacities[1]:
            next_states.append(((a, capacities[1], c), path + ['fill B']))
        if c < capacities[2]:
            next_states.append(((a, b, capacities[2]), path + ['fill C']))
        
        # Empty operations
        if a > 0:
            next_states.append(((0, b, c), path + ['empty A']))
        if b > 0:
            next_states.append(((a, 0, c), path + ['empty B']))
        if c > 0:
            next_states.append(((a, b, 0), path + ['empty C']))
        
        # Pour operations
        # Pour A -> B
        pour_amount = min(a, capacities[1] - b)
        if pour_amount > 0:
            next_states.append(((a - pour_amount, b + pour_amount, c), path + ['pour A->B']))
        
        # Pour A -> C
        pour_amount = min(a, capacities[2] - c)
        if pour_amount > 0:
            next_states.append(((a - pour_amount, b, c + pour_amount), path + ['pour A->C']))
        
        # Pour B -> A
        pour_amount = min(b, capacities[0] - a)
        if pour_amount > 0:
            next_states.append(((a + pour_amount, b - pour_amount, c), path + ['pour B->A']))
        
        # Pour B -> C
        pour_amount = min(b, capacities[2] - c)
        if pour_amount > 0:
            next_states.append(((a, b - pour_amount, c + pour_amount), path + ['pour B->C']))
        
        # Pour C -> A
        pour_amount = min(c, capacities[0] - a)
        if pour_amount > 0:
            next_states.append(((a + pour_amount, b, c - pour_amount), path + ['pour C->A']))
        
        # Pour C -> B
        pour_amount = min(c, capacities[1] - b)
        if pour_amount > 0:
            next_states.append(((a, b + pour_amount, c - pour_amount), path + ['pour C->B']))
        
        # Add new states to the queue
        for state, new_path in next_states:
            if state not in visited:
                visited.add(state)
                queue.append((state, new_path))

# Execute the function and print the result
print(solve_jug_problem())