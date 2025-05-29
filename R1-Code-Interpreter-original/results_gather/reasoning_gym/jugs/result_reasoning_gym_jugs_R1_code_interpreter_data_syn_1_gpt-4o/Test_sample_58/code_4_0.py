from collections import deque

def is_goal(state, target):
    return target in state

def get_next_states(state, capacities):
    next_states = []
    A, B, C = state
    capA, capB, capC = capacities
    
    # Fill operations
    next_states.append(((capA, B, C), "fill A"))  # Fill A
    next_states.append(((A, capB, C), "fill B"))  # Fill B
    next_states.append(((A, B, capC), "fill C"))  # Fill C
    
    # Empty operations
    next_states.append(((0, B, C), "empty A"))  # Empty A
    next_states.append(((A, 0, C), "empty B"))  # Empty B
    next_states.append(((A, B, 0), "empty C"))  # Empty C
    
    # Pour operations
    # Pour A -> B
    pour = min(A, capB - B)
    next_states.append(((A - pour, B + pour, C), "pour A->B"))
    
    # Pour A -> C
    pour = min(A, capC - C)
    next_states.append(((A - pour, B, C + pour), "pour A->C"))
    
    # Pour B -> A
    pour = min(B, capA - A)
    next_states.append(((A + pour, B - pour, C), "pour B->A"))
    
    # Pour B -> C
    pour = min(B, capC - C)
    next_states.append(((A, B - pour, C + pour), "pour B->C"))
    
    # Pour C -> A
    pour = min(C, capA - A)
    next_states.append(((A + pour, B, C - pour), "pour C->A"))
    
    # Pour C -> B
    pour = min(C, capB - B)
    next_states.append(((A, B + pour, C - pour), "pour C->B"))
    
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
        
        for next_state, action in get_next_states(current_state, capacities):
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + [action]))
    
    return None

capacities = (13, 11, 11)
target = 1
solution = bfs(capacities, target)

print(solution)