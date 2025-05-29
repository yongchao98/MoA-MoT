from collections import deque

def is_goal(state):
    return 5 in state

def get_next_states(state):
    A, B, C = state
    next_states = []
    capacities = (10, 13, 13)
    
    # Fill operations
    if A < capacities[0]:
        next_states.append((capacities[0], B, C))
    if B < capacities[1]:
        next_states.append((A, capacities[1], C))
    if C < capacities[2]:
        next_states.append((A, B, capacities[2]))
    
    # Empty operations
    if A > 0:
        next_states.append((0, B, C))
    if B > 0:
        next_states.append((A, 0, C))
    if C > 0:
        next_states.append((A, B, 0))
    
    # Pour operations
    # Pour A -> B
    pour = min(A, capacities[1] - B)
    next_states.append((A - pour, B + pour, C))
    
    # Pour A -> C
    pour = min(A, capacities[2] - C)
    next_states.append((A - pour, B, C + pour))
    
    # Pour B -> A
    pour = min(B, capacities[0] - A)
    next_states.append((A + pour, B - pour, C))
    
    # Pour B -> C
    pour = min(B, capacities[2] - C)
    next_states.append((A, B - pour, C + pour))
    
    # Pour C -> A
    pour = min(C, capacities[0] - A)
    next_states.append((A + pour, B, C - pour))
    
    # Pour C -> B
    pour = min(C, capacities[1] - B)
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
                move = (current_state, next_state)
                queue.append((next_state, path + [move]))
    
    return None

def format_moves(moves):
    formatted_moves = []
    for (A1, B1, C1), (A2, B2, C2) in moves:
        if A1 != A2:
            if A2 > A1:
                formatted_moves.append("fill A")
            else:
                formatted_moves.append("empty A")
        elif B1 != B2:
            if B2 > B1:
                formatted_moves.append("fill B")
            else:
                formatted_moves.append("empty B")
        elif C1 != C2:
            if C2 > C1:
                formatted_moves.append("fill C")
            else:
                formatted_moves.append("empty C")
        else:
            if A1 > A2:
                if B2 > B1:
                    formatted_moves.append("pour A->B")
                else:
                    formatted_moves.append("pour A->C")
            elif B1 > B2:
                if A2 > A1:
                    formatted_moves.append("pour B->A")
                else:
                    formatted_moves.append("pour B->C")
            elif C1 > C2:
                if A2 > A1:
                    formatted_moves.append("pour C->A")
                else:
                    formatted_moves.append("pour C->B")
    
    return formatted_moves

moves = bfs()
formatted_moves = format_moves(moves)
print(formatted_moves)