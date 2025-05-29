from collections import deque

def is_target_reached(state, target):
    return target in state

def get_possible_moves(state, capacities):
    A, B, C = state
    capA, capB, capC = capacities
    moves = []
    
    # Fill operations
    if A < capA:
        moves.append((capA, B, C, 'fill A'))
    if B < capB:
        moves.append((A, capB, C, 'fill B'))
    if C < capC:
        moves.append((A, B, capC, 'fill C'))
    
    # Empty operations
    if A > 0:
        moves.append((0, B, C, 'empty A'))
    if B > 0:
        moves.append((A, 0, C, 'empty B'))
    if C > 0:
        moves.append((A, B, 0, 'empty C'))
    
    # Pour operations
    # Pour A -> B
    if A > 0 and B < capB:
        pour_amount = min(A, capB - B)
        moves.append((A - pour_amount, B + pour_amount, C, 'pour A->B'))
    # Pour A -> C
    if A > 0 and C < capC:
        pour_amount = min(A, capC - C)
        moves.append((A - pour_amount, B, C + pour_amount, 'pour A->C'))
    # Pour B -> A
    if B > 0 and A < capA:
        pour_amount = min(B, capA - A)
        moves.append((A + pour_amount, B - pour_amount, C, 'pour B->A'))
    # Pour B -> C
    if B > 0 and C < capC:
        pour_amount = min(B, capC - C)
        moves.append((A, B - pour_amount, C + pour_amount, 'pour B->C'))
    # Pour C -> A
    if C > 0 and A < capA:
        pour_amount = min(C, capA - A)
        moves.append((A + pour_amount, B, C - pour_amount, 'pour C->A'))
    # Pour C -> B
    if C > 0 and B < capB:
        pour_amount = min(C, capB - B)
        moves.append((A, B + pour_amount, C - pour_amount, 'pour C->B'))
    
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
        
        for next_state, move in get_possible_moves(current_state, capacities):
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + [move]))
    
    return None

capacities = (5, 12, 5)
target = 8
solution = solve_jug_problem(capacities, target)
print(solution)