from collections import deque

def is_target_reached(state, target):
    return target in state

def get_possible_moves(state, capacities):
    A, B, C = state
    maxA, maxB, maxC = capacities
    moves = []
    
    # Fill operations
    if A < maxA:
        moves.append(('fill A', (maxA, B, C)))
    if B < maxB:
        moves.append(('fill B', (A, maxB, C)))
    if C < maxC:
        moves.append(('fill C', (A, B, maxC)))
    
    # Empty operations
    if A > 0:
        moves.append(('empty A', (0, B, C)))
    if B > 0:
        moves.append(('empty B', (A, 0, C)))
    if C > 0:
        moves.append(('empty C', (A, B, 0)))
    
    # Pour operations
    if A > 0 and B < maxB:
        pour_amount = min(A, maxB - B)
        moves.append(('pour A->B', (A - pour_amount, B + pour_amount, C)))
    if A > 0 and C < maxC:
        pour_amount = min(A, maxC - C)
        moves.append(('pour A->C', (A - pour_amount, B, C + pour_amount)))
    if B > 0 and A < maxA:
        pour_amount = min(B, maxA - A)
        moves.append(('pour B->A', (A + pour_amount, B - pour_amount, C)))
    if B > 0 and C < maxC:
        pour_amount = min(B, maxC - C)
        moves.append(('pour B->C', (A, B - pour_amount, C + pour_amount)))
    if C > 0 and A < maxA:
        pour_amount = min(C, maxA - A)
        moves.append(('pour C->A', (A + pour_amount, B, C - pour_amount)))
    if C > 0 and B < maxB:
        pour_amount = min(C, maxB - B)
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

capacities = (4, 13, 4)
target = 7
solution = solve_jug_problem(capacities, target)
print(solution)