from collections import deque

def is_goal(state):
    return 10 in state

def get_possible_moves(state):
    A, B, C = state
    moves = []
    # Fill operations
    if A < 13: moves.append(('fill A', (13, B, C)))
    if B < 6: moves.append(('fill B', (A, 6, C)))
    if C < 13: moves.append(('fill C', (A, B, 13)))
    # Empty operations
    if A > 0: moves.append(('empty A', (0, B, C)))
    if B > 0: moves.append(('empty B', (A, 0, C)))
    if C > 0: moves.append(('empty C', (A, B, 0)))
    # Pour operations
    if A > 0 and B < 6:  # Pour A -> B
        pour_amount = min(A, 6 - B)
        moves.append(('pour A->B', (A - pour_amount, B + pour_amount, C)))
    if A > 0 and C < 13:  # Pour A -> C
        pour_amount = min(A, 13 - C)
        moves.append(('pour A->C', (A - pour_amount, B, C + pour_amount)))
    if B > 0 and A < 13:  # Pour B -> A
        pour_amount = min(B, 13 - A)
        moves.append(('pour B->A', (A + pour_amount, B - pour_amount, C)))
    if B > 0 and C < 13:  # Pour B -> C
        pour_amount = min(B, 13 - C)
        moves.append(('pour B->C', (A, B - pour_amount, C + pour_amount)))
    if C > 0 and A < 13:  # Pour C -> A
        pour_amount = min(C, 13 - A)
        moves.append(('pour C->A', (A + pour_amount, B, C - pour_amount)))
    if C > 0 and B < 6:  # Pour C -> B
        pour_amount = min(C, 6 - B)
        moves.append(('pour C->B', (A, B + pour_amount, C - pour_amount)))
    return moves

def solve_jug_problem():
    initial_state = (0, 0, 0)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, path = queue.popleft()
        if is_goal(current_state):
            return path
        for move, new_state in get_possible_moves(current_state):
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + [move]))

solution = solve_jug_problem()
print(solution)