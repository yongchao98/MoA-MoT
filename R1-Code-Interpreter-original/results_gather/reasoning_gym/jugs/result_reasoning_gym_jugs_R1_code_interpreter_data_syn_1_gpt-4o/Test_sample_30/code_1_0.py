from collections import deque

def is_target_reached(state, target):
    return target in state

def get_next_states(state):
    A, B, C = state
    next_states = []
    # Fill operations
    next_states.append(((5, B, C), 'fill A'))  # Fill A
    next_states.append(((A, 11, C), 'fill B')) # Fill B
    next_states.append(((A, B, 5), 'fill C'))  # Fill C
    # Empty operations
    next_states.append(((0, B, C), 'empty A'))  # Empty A
    next_states.append(((A, 0, C), 'empty B'))  # Empty B
    next_states.append(((A, B, 0), 'empty C'))  # Empty C
    # Pour operations
    # Pour A -> B
    pour = min(A, 11 - B)
    next_states.append(((A - pour, B + pour, C), 'pour A->B'))
    # Pour A -> C
    pour = min(A, 5 - C)
    next_states.append(((A - pour, B, C + pour), 'pour A->C'))
    # Pour B -> A
    pour = min(B, 5 - A)
    next_states.append(((A + pour, B - pour, C), 'pour B->A'))
    # Pour B -> C
    pour = min(B, 5 - C)
    next_states.append(((A, B - pour, C + pour), 'pour B->C'))
    # Pour C -> A
    pour = min(C, 5 - A)
    next_states.append(((A + pour, B, C - pour), 'pour C->A'))
    # Pour C -> B
    pour = min(C, 11 - B)
    next_states.append(((A, B + pour, C - pour), 'pour C->B'))
    return next_states

def get_moves_sequence():
    initial_state = (0, 0, 0)
    target = 3
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, moves = queue.popleft()
        if is_target_reached(current_state, target):
            return moves
        for next_state, move in get_next_states(current_state):
            if next_state not in visited:
                visited.add(next_state)
                new_moves = moves + [move]
                queue.append((next_state, new_moves))
    return []

moves_sequence = get_moves_sequence()
print(moves_sequence)