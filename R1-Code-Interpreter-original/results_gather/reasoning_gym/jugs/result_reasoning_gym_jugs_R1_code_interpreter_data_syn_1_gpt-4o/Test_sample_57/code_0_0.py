from collections import deque

def is_target_reached(state, target):
    return target in state

def get_possible_moves(state, capacities):
    A, B, C = state
    capA, capB, capC = capacities
    moves = []
    
    # Fill operations
    if A < capA:
        moves.append(('fill A', (capA, B, C)))
    if B < capB:
        moves.append(('fill B', (A, capB, C)))
    if C < capC:
        moves.append(('fill C', (A, B, capC)))
    
    # Empty operations
    if A > 0:
        moves.append(('empty A', (0, B, C)))
    if B > 0:
        moves.append(('empty B', (A, 0, C)))
    if C > 0:
        moves.append(('empty C', (A, B, 0)))
    
    # Pour operations
    # Pour A -> B
    if A > 0 and B < capB:
        pour_amount = min(A, capB - B)
        moves.append(('pour A->B', (A - pour_amount, B + pour_amount, C)))
    # Pour A -> C
    if A > 0 and C < capC:
        pour_amount = min(A, capC - C)
        moves.append(('pour A->C', (A - pour_amount, B, C + pour_amount)))
    # Pour B -> A
    if B > 0 and A < capA:
        pour_amount = min(B, capA - A)
        moves.append(('pour B->A', (A + pour_amount, B - pour_amount, C)))
    # Pour B -> C
    if B > 0 and C < capC:
        pour_amount = min(B, capC - C)
        moves.append(('pour B->C', (A, B - pour_amount, C + pour_amount)))
    # Pour C -> A
    if C > 0 and A < capA:
        pour_amount = min(C, capA - A)
        moves.append(('pour C->A', (A + pour_amount, B, C - pour_amount)))
    # Pour C -> B
    if C > 0 and B < capB:
        pour_amount = min(C, capB - B)
        moves.append(('pour C->B', (A, B + pour_amount, C - pour_amount)))
    
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

capacities = (9, 7, 9)
target = 8
solution = solve_jug_problem(capacities, target)
print(solution)