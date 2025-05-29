from collections import deque

def is_target_reached(state, target):
    return target in state

def get_next_states(state):
    A, B, C = state
    next_states = []
    # Fill operations
    next_states.append((6, B, C))  # Fill A
    next_states.append((A, 11, C))  # Fill B
    next_states.append((A, B, 11))  # Fill C
    # Empty operations
    next_states.append((0, B, C))  # Empty A
    next_states.append((A, 0, C))  # Empty B
    next_states.append((A, B, 0))  # Empty C
    # Pour operations
    # Pour A -> B
    pour = min(A, 11 - B)
    next_states.append((A - pour, B + pour, C))
    # Pour A -> C
    pour = min(A, 11 - C)
    next_states.append((A - pour, B, C + pour))
    # Pour B -> A
    pour = min(B, 6 - A)
    next_states.append((A + pour, B - pour, C))
    # Pour B -> C
    pour = min(B, 11 - C)
    next_states.append((A, B - pour, C + pour))
    # Pour C -> A
    pour = min(C, 6 - A)
    next_states.append((A + pour, B, C - pour))
    # Pour C -> B
    pour = min(C, 11 - B)
    next_states.append((A, B + pour, C - pour))
    return next_states

def get_moves(state, next_state):
    A, B, C = state
    A_next, B_next, C_next = next_state
    if A != A_next:
        if A_next == 6:
            return 'fill A'
        elif A_next == 0:
            return 'empty A'
        elif A_next < A:
            return 'pour A->' + ('B' if B_next > B else 'C')
    if B != B_next:
        if B_next == 11:
            return 'fill B'
        elif B_next == 0:
            return 'empty B'
        elif B_next < B:
            return 'pour B->' + ('A' if A_next > A else 'C')
    if C != C_next:
        if C_next == 11:
            return 'fill C'
        elif C_next == 0:
            return 'empty C'
        elif C_next < C:
            return 'pour C->' + ('A' if A_next > A else 'B')
    return ''

def solve_jug_problem(target):
    initial_state = (0, 0, 0)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, moves = queue.popleft()
        if is_target_reached(current_state, target):
            return moves
        for next_state in get_next_states(current_state):
            if next_state not in visited:
                visited.add(next_state)
                move = get_moves(current_state, next_state)
                queue.append((next_state, moves + [move]))

    return []

target = 3
solution = solve_jug_problem(target)
print(solution)