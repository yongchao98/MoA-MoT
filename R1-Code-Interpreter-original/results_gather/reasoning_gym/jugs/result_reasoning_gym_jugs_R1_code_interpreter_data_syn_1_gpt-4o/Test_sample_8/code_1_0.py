from collections import deque

def is_goal(state, target):
    return target in state

def get_next_states(state):
    A, B, C = state
    next_states = []
    moves = []
    # Fill operations
    next_states.append((11, B, C))
    moves.append('fill A')
    next_states.append((A, 6, C))
    moves.append('fill B')
    next_states.append((A, B, 6))
    moves.append('fill C')
    # Empty operations
    next_states.append((0, B, C))
    moves.append('empty A')
    next_states.append((A, 0, C))
    moves.append('empty B')
    next_states.append((A, B, 0))
    moves.append('empty C')
    # Pour operations
    pour = min(A, 6 - B)
    next_states.append((A - pour, B + pour, C))
    moves.append('pour A->B')
    pour = min(A, 6 - C)
    next_states.append((A - pour, B, C + pour))
    moves.append('pour A->C')
    pour = min(B, 11 - A)
    next_states.append((A + pour, B - pour, C))
    moves.append('pour B->A')
    pour = min(B, 6 - C)
    next_states.append((A, B - pour, C + pour))
    moves.append('pour B->C')
    pour = min(C, 11 - A)
    next_states.append((A + pour, B, C - pour))
    moves.append('pour C->A')
    pour = min(C, 6 - B)
    next_states.append((A, B + pour, C - pour))
    moves.append('pour C->B')
    return next_states, moves

def bfs(target):
    initial_state = (0, 0, 0)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal(current_state, target):
            return path
        
        next_states, moves = get_next_states(current_state)
        for next_state, move in zip(next_states, moves):
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + [move]))
    
    return None

# Find the sequence of moves to get 8 liters in one of the jugs
result = bfs(8)
print(result)