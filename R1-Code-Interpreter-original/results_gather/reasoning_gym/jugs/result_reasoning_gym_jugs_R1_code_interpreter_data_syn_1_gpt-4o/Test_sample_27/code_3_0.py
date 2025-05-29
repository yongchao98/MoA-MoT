from collections import deque

def is_goal(state, target):
    return target in state

def get_next_states(state):
    A, B, C = state
    next_states = []
    # Fill operations
    next_states.append(((7, B, C), "fill A"))  # Fill A
    next_states.append(((A, 5, C), "fill B"))  # Fill B
    next_states.append(((A, B, 12), "fill C")) # Fill C
    # Empty operations
    next_states.append(((0, B, C), "empty A"))  # Empty A
    next_states.append(((A, 0, C), "empty B"))  # Empty B
    next_states.append(((A, B, 0), "empty C"))  # Empty C
    # Pour operations
    # Pour A -> B
    pour = min(A, 5 - B)
    next_states.append(((A - pour, B + pour, C), "pour A->B"))
    # Pour A -> C
    pour = min(A, 12 - C)
    next_states.append(((A - pour, B, C + pour), "pour A->C"))
    # Pour B -> A
    pour = min(B, 7 - A)
    next_states.append(((A + pour, B - pour, C), "pour B->A"))
    # Pour B -> C
    pour = min(B, 12 - C)
    next_states.append(((A, B - pour, C + pour), "pour B->C"))
    # Pour C -> A
    pour = min(C, 7 - A)
    next_states.append(((A + pour, B, C - pour), "pour C->A"))
    # Pour C -> B
    pour = min(C, 5 - B)
    next_states.append(((A, B + pour, C - pour), "pour C->B"))
    return next_states

def bfs(target):
    initial_state = (0, 0, 0)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal(current_state, target):
            return path
        
        for next_state, move in get_next_states(current_state):
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + [move]))
    
    return None

moves = bfs(6)
print(moves)