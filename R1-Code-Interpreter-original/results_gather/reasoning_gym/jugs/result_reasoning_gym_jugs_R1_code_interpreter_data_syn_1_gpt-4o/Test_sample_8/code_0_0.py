from collections import deque

def is_goal(state, target):
    return target in state

def get_next_states(state):
    A, B, C = state
    next_states = []
    # Fill operations
    next_states.append((11, B, C))  # Fill A
    next_states.append((A, 6, C))   # Fill B
    next_states.append((A, B, 6))   # Fill C
    # Empty operations
    next_states.append((0, B, C))   # Empty A
    next_states.append((A, 0, C))   # Empty B
    next_states.append((A, B, 0))   # Empty C
    # Pour operations
    # Pour A -> B
    pour = min(A, 6 - B)
    next_states.append((A - pour, B + pour, C))
    # Pour A -> C
    pour = min(A, 6 - C)
    next_states.append((A - pour, B, C + pour))
    # Pour B -> A
    pour = min(B, 11 - A)
    next_states.append((A + pour, B - pour, C))
    # Pour B -> C
    pour = min(B, 6 - C)
    next_states.append((A, B - pour, C + pour))
    # Pour C -> A
    pour = min(C, 11 - A)
    next_states.append((A + pour, B, C - pour))
    # Pour C -> B
    pour = min(C, 6 - B)
    next_states.append((A, B + pour, C - pour))
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
        
        for next_state in get_next_states(current_state):
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + [next_state]))
    
    return None

# Find the sequence of moves to get 8 liters in one of the jugs
result = bfs(8)
print(result)