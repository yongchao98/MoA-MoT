from collections import deque

def solve_jug_problem():
    # Initial state: all jugs are empty
    initial_state = (0, 0, 0)
    # Target amount of water
    target = 10
    # Jug capacities
    capacities = (13, 13, 4)
    
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
        
        # Generate all possible moves
        # Fill operations
        if a < capacities[0]:
            new_state = (capacities[0], b, c)
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + ["fill A"]))
        if b < capacities[1]:
            new_state = (a, capacities[1], c)
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + ["fill B"]))
        if c < capacities[2]:
            new_state = (a, b, capacities[2])
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + ["fill C"]))
        
        # Empty operations
        if a > 0:
            new_state = (0, b, c)
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + ["empty A"]))
        if b > 0:
            new_state = (a, 0, c)
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + ["empty B"]))
        if c > 0:
            new_state = (a, b, 0)
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + ["empty C"]))
        
        # Pour operations
        # Pour A -> B
        if a > 0 and b < capacities[1]:
            pour_amount = min(a, capacities[1] - b)
            new_state = (a - pour_amount, b + pour_amount, c)
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + ["pour A->B"]))
        
        # Pour A -> C
        if a > 0 and c < capacities[2]:
            pour_amount = min(a, capacities[2] - c)
            new_state = (a - pour_amount, b, c + pour_amount)
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + ["pour A->C"]))
        
        # Pour B -> A
        if b > 0 and a < capacities[0]:
            pour_amount = min(b, capacities[0] - a)
            new_state = (a + pour_amount, b - pour_amount, c)
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + ["pour B->A"]))
        
        # Pour B -> C
        if b > 0 and c < capacities[2]:
            pour_amount = min(b, capacities[2] - c)
            new_state = (a, b - pour_amount, c + pour_amount)
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + ["pour B->C"]))
        
        # Pour C -> A
        if c > 0 and a < capacities[0]:
            pour_amount = min(c, capacities[0] - a)
            new_state = (a + pour_amount, b, c - pour_amount)
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + ["pour C->A"]))
        
        # Pour C -> B
        if c > 0 and b < capacities[1]:
            pour_amount = min(c, capacities[1] - b)
            new_state = (a, b + pour_amount, c - pour_amount)
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + ["pour C->B"]))
    
    return []

# Find the solution
solution = solve_jug_problem()
print(solution)