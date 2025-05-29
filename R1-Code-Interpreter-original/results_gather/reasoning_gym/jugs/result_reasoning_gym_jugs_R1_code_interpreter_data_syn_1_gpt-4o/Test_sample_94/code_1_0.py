from collections import deque

def is_target_reached(state, target):
    return target in state

def get_possible_moves(state, capacities):
    A, B, C = state
    capA, capB, capC = capacities
    moves = []
    
    # Fill operations
    if A < capA:
        moves.append(((capA, B, C), 'fill A'))
    if B < capB:
        moves.append(((A, capB, C), 'fill B'))
    if C < capC:
        moves.append(((A, B, capC), 'fill C'))
    
    # Empty operations
    if A > 0:
        moves.append(((0, B, C), 'empty A'))
    if B > 0:
        moves.append(((A, 0, C), 'empty B'))
    if C > 0:
        moves.append(((A, B, 0), 'empty C'))
    
    # Pour operations
    if A > 0 and B < capB:
        pour = min(A, capB - B)
        moves.append(((A - pour, B + pour, C), 'pour A->B'))
    if A > 0 and C < capC:
        pour = min(A, capC - C)
        moves.append(((A - pour, B, C + pour), 'pour A->C'))
    if B > 0 and A < capA:
        pour = min(B, capA - A)
        moves.append(((A + pour, B - pour, C), 'pour B->A'))
    if B > 0 and C < capC:
        pour = min(B, capC - C)
        moves.append(((A, B - pour, C + pour), 'pour B->C'))
    if C > 0 and A < capA:
        pour = min(C, capA - A)
        moves.append(((A + pour, B, C - pour), 'pour C->A'))
    if C > 0 and B < capB:
        pour = min(C, capB - B)
        moves.append(((A, B + pour, C - pour), 'pour C->B'))
    
    return moves

def solve_jug_problem(capacities, target):
    initial_state = (0, 0, 0)
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        state, path = queue.popleft()
        
        if is_target_reached(state, target):
            return path
        
        if state in visited:
            continue
        
        visited.add(state)
        
        for new_state, move in get_possible_moves(state, capacities):
            if new_state not in visited:
                queue.append((new_state, path + [move]))
    
    return None

capacities = (11, 11, 6)
target = 3
solution = solve_jug_problem(capacities, target)
print(solution)