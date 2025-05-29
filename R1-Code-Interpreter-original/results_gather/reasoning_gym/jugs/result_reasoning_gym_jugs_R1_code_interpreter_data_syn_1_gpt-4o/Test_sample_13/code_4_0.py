from collections import deque

def is_goal(state):
    return 12 in state

def get_next_states(state):
    A, B, C = state
    next_states = []
    # Fill operations
    if A < 13:
        next_states.append((13, B, C))
    if B < 11:
        next_states.append((A, 11, C))
    if C < 13:
        next_states.append((A, B, 13))
    # Empty operations
    if A > 0:
        next_states.append((0, B, C))
    if B > 0:
        next_states.append((A, 0, C))
    if C > 0:
        next_states.append((A, B, 0))
    # Pour operations
    # Pour A -> B
    pour = min(A, 11 - B)
    next_states.append((A - pour, B + pour, C))
    # Pour A -> C
    pour = min(A, 13 - C)
    next_states.append((A - pour, B, C + pour))
    # Pour B -> A
    pour = min(B, 13 - A)
    next_states.append((A + pour, B - pour, C))
    # Pour B -> C
    pour = min(B, 13 - C)
    next_states.append((A, B - pour, C + pour))
    # Pour C -> A
    pour = min(C, 13 - A)
    next_states.append((A + pour, B, C - pour))
    # Pour C -> B
    pour = min(C, 11 - B)
    next_states.append((A, B + pour, C - pour))
    
    return next_states

def bfs():
    initial_state = (0, 0, 0)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal(current_state):
            return path
        
        for next_state in get_next_states(current_state):
            if next_state not in visited:
                visited.add(next_state)
                move = determine_move(current_state, next_state)
                queue.append((next_state, path + [move]))

def determine_move(current, next):
    A1, B1, C1 = current
    A2, B2, C2 = next
    if A1 != A2:
        if A2 == 13:
            return 'fill A'
        elif A2 == 0:
            return 'empty A'
        elif A1 > A2:
            if B1 < B2:
                return 'pour A->B'
            else:
                return 'pour A->C'
    elif B1 != B2:
        if B2 == 11:
            return 'fill B'
        elif B2 == 0:
            return 'empty B'
        elif B1 > B2:
            if A1 < A2:
                return 'pour B->A'
            else:
                return 'pour B->C'
    elif C1 != C2:
        if C2 == 13:
            return 'fill C'
        elif C2 == 0:
            return 'empty C'
        elif C1 > C2:
            if A1 < A2:
                return 'pour C->A'
            else:
                return 'pour C->B'

solution = bfs()
print(solution)