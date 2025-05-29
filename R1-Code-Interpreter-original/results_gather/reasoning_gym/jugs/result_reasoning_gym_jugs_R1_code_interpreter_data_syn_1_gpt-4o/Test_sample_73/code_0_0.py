from collections import deque

def is_target_reached(state, target):
    return target in state

def get_possible_moves(state, capacities):
    a, b, c = state
    cap_a, cap_b, cap_c = capacities
    moves = []
    
    # Fill operations
    if a < cap_a:
        moves.append((cap_a, b, c, 'fill A'))
    if b < cap_b:
        moves.append((a, cap_b, c, 'fill B'))
    if c < cap_c:
        moves.append((a, b, cap_c, 'fill C'))
    
    # Empty operations
    if a > 0:
        moves.append((0, b, c, 'empty A'))
    if b > 0:
        moves.append((a, 0, c, 'empty B'))
    if c > 0:
        moves.append((a, b, 0, 'empty C'))
    
    # Pour operations
    if a > 0 and b < cap_b:
        pour_amount = min(a, cap_b - b)
        moves.append((a - pour_amount, b + pour_amount, c, 'pour A->B'))
    if a > 0 and c < cap_c:
        pour_amount = min(a, cap_c - c)
        moves.append((a - pour_amount, b, c + pour_amount, 'pour A->C'))
    if b > 0 and a < cap_a:
        pour_amount = min(b, cap_a - a)
        moves.append((a + pour_amount, b - pour_amount, c, 'pour B->A'))
    if b > 0 and c < cap_c:
        pour_amount = min(b, cap_c - c)
        moves.append((a, b - pour_amount, c + pour_amount, 'pour B->C'))
    if c > 0 and a < cap_a:
        pour_amount = min(c, cap_a - a)
        moves.append((a + pour_amount, b, c - pour_amount, 'pour C->A'))
    if c > 0 and b < cap_b:
        pour_amount = min(c, cap_b - b)
        moves.append((a, b + pour_amount, c - pour_amount, 'pour C->B'))
    
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

capacities = (9, 11, 11)
target = 10
solution = solve_jug_problem(capacities, target)
print(solution)