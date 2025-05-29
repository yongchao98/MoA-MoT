from collections import deque

def is_target_reached(state, target):
    return target in state

def get_possible_moves(state, capacities):
    a, b, c = state
    cap_a, cap_b, cap_c = capacities
    moves = []
    
    # Fill operations
    if a < cap_a:
        moves.append(((cap_a, b, c), 'fill A'))
    if b < cap_b:
        moves.append(((a, cap_b, c), 'fill B'))
    if c < cap_c:
        moves.append(((a, b, cap_c), 'fill C'))
    
    # Empty operations
    if a > 0:
        moves.append(((0, b, c), 'empty A'))
    if b > 0:
        moves.append(((a, 0, c), 'empty B'))
    if c > 0:
        moves.append(((a, b, 0), 'empty C'))
    
    # Pour operations
    # Pour A -> B
    pour_ab = min(a, cap_b - b)
    if pour_ab > 0:
        moves.append(((a - pour_ab, b + pour_ab, c), 'pour A->B'))
    
    # Pour A -> C
    pour_ac = min(a, cap_c - c)
    if pour_ac > 0:
        moves.append(((a - pour_ac, b, c + pour_ac), 'pour A->C'))
    
    # Pour B -> A
    pour_ba = min(b, cap_a - a)
    if pour_ba > 0:
        moves.append(((a + pour_ba, b - pour_ba, c), 'pour B->A'))
    
    # Pour B -> C
    pour_bc = min(b, cap_c - c)
    if pour_bc > 0:
        moves.append(((a, b - pour_bc, c + pour_bc), 'pour B->C'))
    
    # Pour C -> A
    pour_ca = min(c, cap_a - a)
    if pour_ca > 0:
        moves.append(((a + pour_ca, b, c - pour_ca), 'pour C->A'))
    
    # Pour C -> B
    pour_cb = min(c, cap_b - b)
    if pour_cb > 0:
        moves.append(((a, b + pour_cb, c - pour_cb), 'pour C->B'))
    
    return moves

def find_solution(capacities, target):
    initial_state = (0, 0, 0)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_target_reached(current_state, target):
            return path
        
        for new_state, move in get_possible_moves(current_state, capacities):
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + [move]))
    
    return None

capacities = (7, 9, 9)
target = 6
solution = find_solution(capacities, target)
print(solution)