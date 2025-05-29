from collections import deque

def is_goal(state, target):
    return target in state

def get_next_states(state, capacities):
    A, B, C = state
    capA, capB, capC = capacities
    next_states = []
    
    # Fill operations
    if A < capA:
        next_states.append(((capA, B, C), 'fill A'))
    if B < capB:
        next_states.append(((A, capB, C), 'fill B'))
    if C < capC:
        next_states.append(((A, B, capC), 'fill C'))
    
    # Empty operations
    if A > 0:
        next_states.append(((0, B, C), 'empty A'))
    if B > 0:
        next_states.append(((A, 0, C), 'empty B'))
    if C > 0:
        next_states.append(((A, B, 0), 'empty C'))
    
    # Pour operations
    # Pour A -> B
    pour = min(A, capB - B)
    next_states.append(((A - pour, B + pour, C), 'pour A->B'))
    
    # Pour A -> C
    pour = min(A, capC - C)
    next_states.append(((A - pour, B, C + pour), 'pour A->C'))
    
    # Pour B -> A
    pour = min(B, capA - A)
    next_states.append(((A + pour, B - pour, C), 'pour B->A'))
    
    # Pour B -> C
    pour = min(B, capC - C)
    next_states.append(((A, B - pour, C + pour), 'pour B->C'))
    
    # Pour C -> A
    pour = min(C, capA - A)
    next_states.append(((A + pour, B, C - pour), 'pour C->A'))
    
    # Pour C -> B
    pour = min(C, capB - B)
    next_states.append(((A, B + pour, C - pour), 'pour C->B'))
    
    return next_states

def bfs(capacities, target):
    initial_state = (0, 0, 0)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal(current_state, target):
            return path
        
        for next_state, move in get_next_states(current_state, capacities):
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + [move]))
    
    return None

capacities = (13, 13, 4)
target = 6
solution = bfs(capacities, target)

print(solution)