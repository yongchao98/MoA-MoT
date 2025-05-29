from collections import deque

def is_goal(state):
    return 4 in state

def get_next_states(state):
    A, B, C = state
    next_states = []
    # Fill operations
    next_states.append((13, B, C))  # Fill A
    next_states.append((A, 7, C))   # Fill B
    next_states.append((A, B, 13))  # Fill C
    # Empty operations
    next_states.append((0, B, C))   # Empty A
    next_states.append((A, 0, C))   # Empty B
    next_states.append((A, B, 0))   # Empty C
    # Pour operations
    # Pour A -> B
    pour = min(A, 7 - B)
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
    pour = min(C, 7 - B)
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
    for (from_state, to_state) in moves:
        A1, B1, C1 = from_state
        A2, B2, C2 = to_state
        if A1 != A2:
            if A2 == 13:
                formatted_moves.append("fill A")
            elif A2 == 0:
                formatted_moves.append("empty A")
            else:
                if B1 != B2:
                    formatted_moves.append("pour B->A")
                else:
                    formatted_moves.append("pour C->A")
        elif B1 != B2:
            if B2 == 7:
                formatted_moves.append("fill B")
            elif B2 == 0:
                formatted_moves.append("empty B")
            else:
                if A1 != A2:
                    formatted_moves.append("pour A->B")
                else:
                    formatted_moves.append("pour C->B")
        elif C1 != C2:
            if C2 == 13:
                formatted_moves.append("fill C")
            elif C2 == 0:
                formatted_moves.append("empty C")
            else:
                if A1 != A2:
                    formatted_moves.append("pour A->C")
                else:
                    formatted_moves.append("pour B->C")
    return formatted_moves

moves = bfs()
formatted_moves = format_moves(moves)
print(formatted_moves)