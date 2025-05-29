from collections import deque

def solve_jug_problem():
    # Initial state: all jugs are empty
    initial_state = (0, 0, 0)
    # Target amount of water
    target = 12
    # Jug capacities
    capacities = (11, 13, 11)
    
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
        if a < capacities[0]:
            next_states.append(((capacities[0], b, c), moves + ['fill A']))
        if b < capacities[1]:
            next_states.append(((a, capacities[1], c), moves + ['fill B']))
        if c < capacities[2]:
            next_states.append(((a, b, capacities[2]), moves + ['fill C']))
        
        # Empty operations
        if a > 0:
            next_states.append(((0, b, c), moves + ['empty A']))
        if b > 0:
            next_states.append(((a, 0, c), moves + ['empty B']))
        if c > 0:
            next_states.append(((a, b, 0), moves + ['empty C']))
        
        # Pour operations
        # Pour A -> B
        if a > 0 and b < capacities[1]:
            pour_amount = min(a, capacities[1] - b)
            next_states.append(((a - pour_amount, b + pour_amount, c), moves + ['pour A->B']))
        
        # Pour A -> C
        if a > 0 and c < capacities[2]:
            pour_amount = min(a, capacities[2] - c)
            next_states.append(((a - pour_amount, b, c + pour_amount), moves + ['pour A->C']))
        
        # Pour B -> A
        if b > 0 and a < capacities[0]:
            pour_amount = min(b, capacities[0] - a)
            next_states.append(((a + pour_amount, b - pour_amount, c), moves + ['pour B->A']))
        
        # Pour B -> C
        if b > 0 and c < capacities[2]:
            pour_amount = min(b, capacities[2] - c)
            next_states.append(((a, b - pour_amount, c + pour_amount), moves + ['pour B->C']))
        
        # Pour C -> A
        if c > 0 and a < capacities[0]:
            pour_amount = min(c, capacities[0] - a)
            next_states.append(((a + pour_amount, b, c - pour_amount), moves + ['pour C->A']))
        
        # Pour C -> B
        if c > 0 and b < capacities[1]:
            pour_amount = min(c, capacities[1] - b)
            next_states.append(((a, b + pour_amount, c - pour_amount), moves + ['pour C->B']))
        
        # Add valid next states to the queue
        for state, new_moves in next_states:
            if state not in visited:
                visited.add(state)
                queue.append((state, new_moves))

# Execute the function and print the result
result = solve_jug_problem()
print(result)