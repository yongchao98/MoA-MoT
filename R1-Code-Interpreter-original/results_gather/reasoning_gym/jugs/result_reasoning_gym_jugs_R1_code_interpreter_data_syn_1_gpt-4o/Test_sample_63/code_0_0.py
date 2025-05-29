from collections import deque

def is_goal(state, target):
    return target in state

def get_next_states(state, capacities):
    A, B, C = state
    capA, capB, capC = capacities
    next_states = []
    
    # Fill operations
    if A < capA:
        next_states.append((capA, B, C))
    if B < capB:
        next_states.append((A, capB, C))
    if C < capC:
        next_states.append((A, B, capC))
    
    # Empty operations
    if A > 0:
        next_states.append((0, B, C))
    if B > 0:
        next_states.append((A, 0, C))
    if C > 0:
        next_states.append((A, B, 0))
    
    # Pour operations
    # Pour A -> B
    pour = min(A, capB - B)
    next_states.append((A - pour, B + pour, C))
    
    # Pour A -> C
    pour = min(A, capC - C)
    next_states.append((A - pour, B, C + pour))
    
    # Pour B -> A
    pour = min(B, capA - A)
    next_states.append((A + pour, B - pour, C))
    
    # Pour B -> C
    pour = min(B, capC - C)
    next_states.append((A, B - pour, C + pour))
    
    # Pour C -> A
    pour = min(C, capA - A)
    next_states.append((A + pour, B, C - pour))
    
    # Pour C -> B
    pour = min(C, capB - B)
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

capacities = (13, 13, 4)
target = 6
solution = bfs(capacities, target)

# Convert the solution to a JSON-parsable list of moves
moves = []
for (start, end) in solution:
    if start[0] != end[0]:
        if end[0] == capacities[0]:
            moves.append("fill A")
        elif end[0] == 0:
            moves.append("empty A")
        else:
            if start[1] != end[1]:
                moves.append("pour B->A")
            else:
                moves.append("pour C->A")
    elif start[1] != end[1]:
        if end[1] == capacities[1]:
            moves.append("fill B")
        elif end[1] == 0:
            moves.append("empty B")
        else:
            if start[0] != end[0]:
                moves.append("pour A->B")
            else:
                moves.append("pour C->B")
    elif start[2] != end[2]:
        if end[2] == capacities[2]:
            moves.append("fill C")
        elif end[2] == 0:
            moves.append("empty C")
        else:
            if start[0] != end[0]:
                moves.append("pour A->C")
            else:
                moves.append("pour B->C")

print(moves)