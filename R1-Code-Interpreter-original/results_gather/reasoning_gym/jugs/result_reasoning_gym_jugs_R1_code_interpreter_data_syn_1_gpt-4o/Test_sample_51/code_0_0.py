from collections import deque

def solve_jug_problem():
    # Initial state: all jugs are empty
    initial_state = (0, 0, 0)
    # Target amount of water
    target = 7
    # Jug capacities
    capacities = (8, 13, 13)
    
    # Queue for BFS
    queue = deque([(initial_state, [])])
    # Set to keep track of visited states
    visited = set()
    visited.add(initial_state)
    
    while queue:
        (a, b, c), moves = queue.popleft()
        
        # Check if we have reached the target
        if a == target or b == target or c == target:
            return moves
        
        # Generate all possible next states
        next_states = []
        
        # Fill operations
        next_states.append(((capacities[0], b, c), moves + ['fill A']))
        next_states.append(((a, capacities[1], c), moves + ['fill B']))
        next_states.append(((a, b, capacities[2]), moves + ['fill C']))
        
        # Empty operations
        next_states.append(((0, b, c), moves + ['empty A']))
        next_states.append(((a, 0, c), moves + ['empty B']))
        next_states.append(((a, b, 0), moves + ['empty C']))
        
        # Pour operations
        # Pour A -> B
        pour = min(a, capacities[1] - b)
        next_states.append(((a - pour, b + pour, c), moves + ['pour A->B']))
        
        # Pour A -> C
        pour = min(a, capacities[2] - c)
        next_states.append(((a - pour, b, c + pour), moves + ['pour A->C']))
        
        # Pour B -> A
        pour = min(b, capacities[0] - a)
        next_states.append(((a + pour, b - pour, c), moves + ['pour B->A']))
        
        # Pour B -> C
        pour = min(b, capacities[2] - c)
        next_states.append(((a, b - pour, c + pour), moves + ['pour B->C']))
        
        # Pour C -> A
        pour = min(c, capacities[0] - a)
        next_states.append(((a + pour, b, c - pour), moves + ['pour C->A']))
        
        # Pour C -> B
        pour = min(c, capacities[1] - b)
        next_states.append(((a, b + pour, c - pour), moves + ['pour C->B']))
        
        # Add valid next states to the queue
        for state, new_moves in next_states:
            if state not in visited:
                visited.add(state)
                queue.append((state, new_moves))

# Execute the function and print the result
print(solve_jug_problem())