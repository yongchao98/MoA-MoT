from collections import deque

def is_target_reached(state, target):
    return target in state

def get_possible_moves(state, capacities):
    moves = []
    a, b, c = state
    cap_a, cap_b, cap_c = capacities
    
    # Fill operations
    if a < cap_a:
        moves.append(('fill A', (cap_a, b, c)))
    if b < cap_b:
        moves.append(('fill B', (a, cap_b, c)))
    if c < cap_c:
        moves.append(('fill C', (a, b, cap_c)))
    
    # Empty operations
    if a > 0:
        moves.append(('empty A', (0, b, c)))
    if b > 0:
        moves.append(('empty B', (a, 0, c)))
    if c > 0:
        moves.append(('empty C', (a, b, 0)))
    
    # Pour operations
    # Pour A -> B
    if a > 0 and b < cap_b:
        pour_amount = min(a, cap_b - b)
        moves.append(('pour A->B', (a - pour_amount, b + pour_amount, c)))
    # Pour A -> C
    if a > 0 and c < cap_c:
        pour_amount = min(a, cap_c - c)
        moves.append(('pour A->C', (a - pour_amount, b, c + pour_amount)))
    # Pour B -> A
    if b > 0 and a < cap_a:
        pour_amount = min(b, cap_a - a)
        moves.append(('pour B->A', (a + pour_amount, b - pour_amount, c)))
    # Pour B -> C
    if b > 0 and c < cap_c:
        pour_amount = min(b, cap_c - c)
        moves.append(('pour B->C', (a, b - pour_amount, c + pour_amount)))
    # Pour C -> A
    if c > 0 and a < cap_a:
        pour_amount = min(c, cap_a - a)
        moves.append(('pour C->A', (a + pour_amount, b, c - pour_amount)))
    # Pour C -> B
    if c > 0 and b < cap_b:
        pour_amount = min(c, cap_b - b)
        moves.append(('pour C->B', (a, b + pour_amount, c - pour_amount)))
    
    return moves

def solve_jug_problem(capacities, target):
    initial_state = (0, 0, 0)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_target_reached(current_state, target):
            return path
        
        for move, new_state in get_possible_moves(current_state, capacities):
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + [move]))
    
    return None

capacities = (6, 7, 7)
target = 3
solution = solve_jug_problem(capacities, target)
print(solution)