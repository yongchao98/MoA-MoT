from collections import deque

def get_possible_moves(state):
    A, B, C = state
    moves = []
    # Fill operations
    if A < 4:
        moves.append(('fill A', (4, B, C)))
    if B < 4:
        moves.append(('fill B', (A, 4, C)))
    if C < 9:
        moves.append(('fill C', (A, B, 9)))
    
    # Empty operations
    if A > 0:
        moves.append(('empty A', (0, B, C)))
    if B > 0:
        moves.append(('empty B', (A, 0, C)))
    if C > 0:
        moves.append(('empty C', (A, B, 0)))
    
    # Pour operations
    if A > 0 and B < 4:
        pour_amount = min(A, 4 - B)
        moves.append(('pour A->B', (A - pour_amount, B + pour_amount, C)))
    if A > 0 and C < 9:
        pour_amount = min(A, 9 - C)
        moves.append(('pour A->C', (A - pour_amount, B, C + pour_amount)))
    if B > 0 and A < 4:
        pour_amount = min(B, 4 - A)
        moves.append(('pour B->A', (A + pour_amount, B - pour_amount, C)))
    if B > 0 and C < 9:
        pour_amount = min(B, 9 - C)
        moves.append(('pour B->C', (A, B - pour_amount, C + pour_amount)))
    if C > 0 and A < 4:
        pour_amount = min(C, 4 - A)
        moves.append(('pour C->A', (A + pour_amount, B, C - pour_amount)))
    if C > 0 and B < 4:
        pour_amount = min(C, 4 - B)
        moves.append(('pour C->B', (A, B + pour_amount, C - pour_amount)))
    
    return moves

def find_solution():
    initial_state = (0, 0, 0)
    target = 7
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        if target in current_state:
            return path
        
        for move, new_state in get_possible_moves(current_state):
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + [move]))
    
    return None

solution = find_solution()
print(solution)