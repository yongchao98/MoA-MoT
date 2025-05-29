from collections import deque

def is_target_reached(state, target):
    return target in state

def get_possible_moves(state, capacities):
    A, B, C = state
    maxA, maxB, maxC = capacities
    moves = []
    
    # Fill operations
    if A < maxA:
        moves.append(((maxA, B, C), 'fill A'))
    if B < maxB:
        moves.append(((A, maxB, C), 'fill B'))
    if C < maxC:
        moves.append(((A, B, maxC), 'fill C'))
    
    # Empty operations
    if A > 0:
        moves.append(((0, B, C), 'empty A'))
    if B > 0:
        moves.append(((A, 0, C), 'empty B'))
    if C > 0:
        moves.append(((A, B, 0), 'empty C'))
    
    # Pour operations
    # Pour A -> B
    if A > 0 and B < maxB:
        pour_amount = min(A, maxB - B)
        moves.append(((A - pour_amount, B + pour_amount, C), 'pour A->B'))
    # Pour A -> C
    if A > 0 and C < maxC:
        pour_amount = min(A, maxC - C)
        moves.append(((A - pour_amount, B, C + pour_amount), 'pour A->C'))
    # Pour B -> A
    if B > 0 and A < maxA:
        pour_amount = min(B, maxA - A)
        moves.append(((A + pour_amount, B - pour_amount, C), 'pour B->A'))
    # Pour B -> C
    if B > 0 and C < maxC:
        pour_amount = min(B, maxC - C)
        moves.append(((A, B - pour_amount, C + pour_amount), 'pour B->C'))
    # Pour C -> A
    if C > 0 and A < maxA:
        pour_amount = min(C, maxA - A)
        moves.append(((A + pour_amount, B, C - pour_amount), 'pour C->A'))
    # Pour C -> B
    if C > 0 and B < maxB:
        pour_amount = min(C, maxB - B)
        moves.append(((A, B + pour_amount, C - pour_amount), 'pour C->B'))
    
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

capacities = (10, 10, 9)
target = 4
solution = solve_jug_problem(capacities, target)
print(solution)