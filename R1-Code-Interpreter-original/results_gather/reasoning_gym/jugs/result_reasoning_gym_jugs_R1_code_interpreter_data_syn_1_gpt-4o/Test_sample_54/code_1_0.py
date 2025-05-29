from collections import deque

def is_goal(state, target):
    return target in state

def get_next_states(state, capacities):
    A, B, C = state
    maxA, maxB, maxC = capacities
    next_states = []
    
    # Fill operations
    next_states.append((maxA, B, C))  # Fill A
    next_states.append((A, maxB, C))  # Fill B
    next_states.append((A, B, maxC))  # Fill C
    
    # Empty operations
    next_states.append((0, B, C))  # Empty A
    next_states.append((A, 0, C))  # Empty B
    next_states.append((A, B, 0))  # Empty C
    
    # Pour operations
    # Pour A -> B
    pour = min(A, maxB - B)
    next_states.append((A - pour, B + pour, C))
    
    # Pour A -> C
    pour = min(A, maxC - C)
    next_states.append((A - pour, B, C + pour))
    
    # Pour B -> A
    pour = min(B, maxA - A)
    next_states.append((A + pour, B - pour, C))
    
    # Pour B -> C
    pour = min(B, maxC - C)
    next_states.append((A, B - pour, C + pour))
    
    # Pour C -> A
    pour = min(C, maxA - A)
    next_states.append((A + pour, B, C - pour))
    
    # Pour C -> B
    pour = min(C, maxB - B)
    next_states.append((A, B + pour, C - pour))
    
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
        
        for next_state in get_next_states(current_state, capacities):
            if next_state not in visited:
                visited.add(next_state)
                move = (current_state, next_state)
                queue.append((next_state, path + [move]))
    
    return None

capacities = (11, 13, 9)
target = 12
solution = bfs(capacities, target)

# Convert the solution to a list of moves
moves = []
for (from_state, to_state) in solution:
    A1, B1, C1 = from_state
    A2, B2, C2 = to_state
    if A1 != A2:
        if A2 == capacities[0]:
            moves.append("fill A")
        elif A2 == 0:
            moves.append("empty A")
        else:
            if A1 > A2:
                moves.append("pour A->B" if B1 != B2 else "pour A->C")
            else:
                moves.append("pour B->A" if B1 != B2 else "pour C->A")
    elif B1 != B2:
        if B2 == capacities[1]:
            moves.append("fill B")
        elif B2 == 0:
            moves.append("empty B")
        else:
            if B1 > B2:
                moves.append("pour B->A" if A1 != A2 else "pour B->C")
            else:
                moves.append("pour A->B" if A1 != A2 else "pour C->B")
    elif C1 != C2:
        if C2 == capacities[2]:
            moves.append("fill C")
        elif C2 == 0:
            moves.append("empty C")
        else:
            if C1 > C2:
                moves.append("pour C->A" if A1 != A2 else "pour C->B")
            else:
                moves.append("pour A->C" if A1 != A2 else "pour B->C")

print(moves)