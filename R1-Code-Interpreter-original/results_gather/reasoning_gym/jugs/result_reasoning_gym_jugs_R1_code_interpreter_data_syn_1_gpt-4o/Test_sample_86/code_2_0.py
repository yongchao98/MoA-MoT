from collections import deque

def is_goal(state, target):
    return target in state

def get_next_states(state):
    A, B, C = state
    next_states = []
    # Fill operations
    next_states.append((13, B, C))  # Fill A
    next_states.append((A, 3, C))   # Fill B
    next_states.append((A, B, 3))   # Fill C
    # Empty operations
    next_states.append((0, B, C))   # Empty A
    next_states.append((A, 0, C))   # Empty B
    next_states.append((A, B, 0))   # Empty C
    # Pour operations
    # Pour A -> B
    pour = min(A, 3 - B)
    next_states.append((A - pour, B + pour, C))
    # Pour A -> C
    pour = min(A, 3 - C)
    next_states.append((A - pour, B, C + pour))
    # Pour B -> A
    pour = min(B, 13 - A)
    next_states.append((A + pour, B - pour, C))
    # Pour B -> C
    pour = min(B, 3 - C)
    next_states.append((A, B - pour, C + pour))
    # Pour C -> A
    pour = min(C, 13 - A)
    next_states.append((A + pour, B, C - pour))
    # Pour C -> B
    pour = min(C, 3 - B)
    next_states.append((A, B + pour, C - pour))
    return next_states

def bfs(initial_state, target):
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
                move = (current_state, next_state)
                queue.append((next_state, path + [move]))
    
    return None

initial_state = (0, 0, 0)
target = 11
solution = bfs(initial_state, target)

# Convert the solution to a list of moves
moves = []
for (from_state, to_state) in solution:
    A1, B1, C1 = from_state
    A2, B2, C2 = to_state
    if A1 != A2:
        if A2 == 13:
            moves.append("fill A")
        elif A2 == 0:
            moves.append("empty A")
        else:
            if B1 != B2:
                moves.append("pour B->A")
            else:
                moves.append("pour C->A")
    elif B1 != B2:
        if B2 == 3:
            moves.append("fill B")
        elif B2 == 0:
            moves.append("empty B")
        else:
            if A1 != A2:
                moves.append("pour A->B")
            else:
                moves.append("pour C->B")
    elif C1 != C2:
        if C2 == 3:
            moves.append("fill C")
        elif C2 == 0:
            moves.append("empty C")
        else:
            if A1 != A2:
                moves.append("pour A->C")
            else:
                moves.append("pour B->C")

print(moves)