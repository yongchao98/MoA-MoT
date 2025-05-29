from collections import deque

def solve_jug_problem():
    # Initial state: all jugs are empty
    initial_state = (0, 0, 0)
    target = 8
    capacities = (7, 7, 10)
    
    # Queue for BFS: each element is (state, moves)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        (a, b, c), moves = queue.popleft()
        
        # Check if any jug has the target amount
        if a == target or b == target or c == target:
            return moves
        
        # Generate all possible new states
        new_states = []
        
        # Fill operations
        new_states.append(((capacities[0], b, c), moves + ['fill A']))
        new_states.append(((a, capacities[1], c), moves + ['fill B']))
        new_states.append(((a, b, capacities[2]), moves + ['fill C']))
        
        # Empty operations
        new_states.append(((0, b, c), moves + ['empty A']))
        new_states.append(((a, 0, c), moves + ['empty B']))
        new_states.append(((a, b, 0), moves + ['empty C']))
        
        # Pour operations
        # Pour A -> B
        pour = min(a, capacities[1] - b)
        new_states.append(((a - pour, b + pour, c), moves + ['pour A->B']))
        
        # Pour A -> C
        pour = min(a, capacities[2] - c)
        new_states.append(((a - pour, b, c + pour), moves + ['pour A->C']))
        
        # Pour B -> A
        pour = min(b, capacities[0] - a)
        new_states.append(((a + pour, b - pour, c), moves + ['pour B->A']))
        
        # Pour B -> C
        pour = min(b, capacities[2] - c)
        new_states.append(((a, b - pour, c + pour), moves + ['pour B->C']))
        
        # Pour C -> A
        pour = min(c, capacities[0] - a)
        new_states.append(((a + pour, b, c - pour), moves + ['pour C->A']))
        
        # Pour C -> B
        pour = min(c, capacities[1] - b)
        new_states.append(((a, b + pour, c - pour), moves + ['pour C->B']))
        
        # Add new states to the queue if not visited
        for new_state, new_moves in new_states:
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, new_moves))

# Execute the function and print the result
result = solve_jug_problem()
print(result)